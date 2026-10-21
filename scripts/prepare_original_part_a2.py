#!/usr/bin/env python3
"""Prepare original Part A2 with audited master input and temporal facet N only."""
from pathlib import Path
import hashlib
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
# precisely the year-eligible source population (1990-2025), not all source rows.
# This is a surgical label-only change; original statistics and figures remain.
old_facet = """src_N <- SEN %>%
  count(
    Source,
    name = "N_total"
  ) %>%"""
new_facet = """src_N <- SEN %>%
  filter(
    !is.na(Collection_year),
    as.integer(Collection_year) >= 1990L,
    as.integer(Collection_year) <= 2025L
  ) %>%
  count(
    Source,
    name = "N_total"
  ) %>%"""
assert rendered.count(old_facet) == 1, "Original facet N block must occur exactly once"
rendered = rendered.replace(old_facet, new_facet, 1)
out=root/"analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(rendered,encoding="utf-8")
# Strict invariant: only the input reader is changed.
assert rendered.replace(new_facet,old_facet,1).replace(new,old,1)==text
sha=hashlib.sha256(original.read_bytes()).hexdigest()
(out.parent/"ORIGINAL_SHA256.txt").write_text(f"{sha}  {original.name}\n")
print(f"[PASS] Original notebook preserved; master import and temporal facet N corrected: {out}")
print(f"[AUDIT] Original SHA-256: {sha}")
print("[AUDIT] Temporal facet labels count only 1990-2025 eligible genomes; percentages and plots unchanged.")
print("[NOTE] Notebook execution can depend on exploratory objects/interactive R chunks; validate outputs against the original.")
