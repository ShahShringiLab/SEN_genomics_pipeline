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
# The original Top-100 chart is horizontal (coord_flip): scale_y is its x-axis.
chart_axis_old = "limits = c(0, 22000),"
chart_axis_new = "limits = c(0, 25000),"
assert rendered.count(chart_axis_old) == 1, "FQR Top-100 axis block must be unique"
rendered = rendered.replace(chart_axis_old, chart_axis_new, 1)

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
# Targeted graphical change: after coord_flip(), the count scale is the
# displayed horizontal (x) axis in FQR_marker_counts_top100.
# Retain the source notebook unchanged; modify only the generated working Rmd.
def set_top100_xmax(doc):
    start = "p_tokens <- ggplot(plot_df, aes(x = token_pretty, y = n_isolates, fill = Type)) +"
    end = "print(p_tokens)"
    assert doc.count(start) == 1, "Expected one Top100 marker chart definition"
    prefix, rest = doc.split(start, 1)
    assert rest.count(end) >= 1, "Cannot find end of Top100 chart"
    plot, suffix = rest.split(end, 1)
    old_limit = "limits = c(0, 22000)"
    assert plot.count(old_limit) == 1, "Expected one Top100 count axis with max 22000"
    plot = plot.replace(old_limit, "limits = c(0, 25000)", 1)
    return prefix + start + plot + end + suffix

rendered = set_top100_xmax(rendered)
# Strict reproducibility: only input, facet N, temporal cutoff and Top100 axis may change.
expected = text.replace(old, new, 1).replace(old_facet, new_facet, 1).replace(chart_axis_old, chart_axis_new, 1)
expected, expected_n = re.subn(lower_cutoff, "", expected)
expected = set_top100_xmax(expected)
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
print("[AUDIT] Full-master Top-100 FQR chart x-axis: 0–25,000.")
print("[AUDIT] FQR_marker_counts_top100 x-axis limit set to 25,000 (after coord_flip).")
print("[NOTE] Notebook execution can depend on exploratory objects/interactive R chunks; validate outputs against the original.")
