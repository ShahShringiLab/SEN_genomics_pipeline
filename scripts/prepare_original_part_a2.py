#!/usr/bin/env python3
"""Prepare original Part A2 with audited master input and consistent temporal eligibility."""
from pathlib import Path
import hashlib
import re
import sys
root=Path(__file__).resolve().parents[1]
original=Path(sys.argv[1]) if len(sys.argv)>1 else root/"provenance/SEN_Genomes2_uploaded_original_20261021.Rmd"
if not original.is_file(): raise SystemExit(f"Missing original Rmd: {original}")
text=original.read_text(encoding="utf-8-sig")
old='SEN <- read_csv(file_path, col_types = cols(.default = col_character()))'
assert text.count(old)==1, f"Expected exactly one original import, saw {text.count(old)}"
# Stage 0 parses and saves 33,831 x 51 before any derived columns. The notebook
# reads ONLY that snapshot; original analysis chunks remain verbatim.
new='''file_path <- file.path(Sys.getenv("SEN_ROOT", unset=getwd()),
                       "results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv")
SEN <- read_csv(file_path, col_types = cols(.default = col_character()))
stopifnot(nrow(SEN)==33831L,ncol(SEN)==51L)'''
rendered=text.replace(old,new,1)
# Original phenotype script belongs to the same notebook. Correct the input
# path ONLY; the original statistical analysis code remains unchanged.
phenotype_old = '"FDASEN_AMRPhenotype2.csv"'
phenotype_new = 'file.path(Sys.getenv("SEN_ROOT", unset=getwd()), "data/part_a/phenotype/FDASEN_AMRPhenotype2.csv")'
assert rendered.count(phenotype_old) == 1, "Original phenotype input path not uniquely identified"
rendered = rendered.replace(phenotype_old, phenotype_new, 1)
# The original dot plot uses source-specific facet labels. Its N must reflect
# precisely the plotted source population: any valid year <=2025.
# Historical years <=2000 are grouped under the single <=2000 category.
old_facet = """src_N <- SEN %>%
  count(
    Source,
    name = "N_total"
  ) %>%"""
new_facet = """src_N <- SEN %>%
  filter(
    !is.na(Collection_year),
    as.integer(Collection_year) <= 2025L
  ) %>%
  count(
    Source,
    name = "N_total"
  ) %>%"""
assert rendered.count(old_facet) == 1, "Original facet N block must occur exactly once"
rendered = rendered.replace(old_facet, new_facet, 1)

# The Top-100 figure uses coord_flip(). Extend the horizontal count axis
# to 250,000; suppress the endpoint's numeric tick label by stopping labeled
# breaks at 200,000. This is an intentional tenfold change from 25,000.
old_top100 = """  scale_y_continuous(
    limits = c(0, 22000),
    expand = expansion(mult = c(0, 0.02))
  ) +"""
new_top100 = """  scale_y_continuous(
    limits = c(0, 250000),
    breaks = seq(0, 200000, by = 50000),
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0.02))
  ) +"""
assert rendered.count(old_top100) == 1, (
    f"Expected one FQR_marker_counts_top100 axis block; found {rendered.count(old_top100)}"
)
rendered = rendered.replace(old_top100, new_top100, 1)

# All six temporal uses (uploaded 30-chunk original) or four uses (archived
# 26-chunk notebook) incorrectly discarded dated genomes before 1990, although
# both notebooks have a <=2000 bin. Remove only those lower-year cutoff LINES
# from the generated working copy. Keep the archived/uploaded original pristine.
lower_cutoff = r"(?m)^[ \t]*Collection_year[ \t]*>=[ \t]*1990[L]?,[ \t]*\n"
rendered, n_dropped = re.subn(lower_cutoff, "", rendered)
assert n_dropped in (4, 6), (
    f"Expected 4 or 6 temporal lower-year cutoffs, found {n_dropped}. "
    "Refusing an unverified edit."
)
# Strict reproducibility: only both input paths, facet N, temporal cutoff and Top100 axis may change.
expected = text.replace(old, new, 1).replace(phenotype_old, phenotype_new, 1).replace(old_facet, new_facet, 1).replace(old_top100, new_top100, 1)
expected, expected_n = re.subn(lower_cutoff, "", expected)
assert expected_n == n_dropped and rendered == expected
out=root/"analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(rendered,encoding="utf-8")
sha=hashlib.sha256(original.read_bytes()).hexdigest()
(out.parent/"ORIGINAL_SHA256.txt").write_text(f"{sha}  {original.name}\n")
print(f"[PASS] Original notebook preserved; master import and temporal facet N corrected: {out}")
print(f"[AUDIT] Original SHA-256: {sha}")
print(f"[AUDIT] Removed {n_dropped} unintended pre-1990 temporal cutoffs; years <=2000 remain grouped.")
print("[AUDIT] Dot plot facets and numerators use the same valid-year <=2025 population.")
print("[AUDIT] FQR_marker_counts_top100 horizontal count axis extends to 250,000; endpoint 250,000 has no tick label.")
print("[NOTE] Notebook execution can depend on exploratory objects/interactive R chunks; validate outputs against the original.")
