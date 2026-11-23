#!/usr/bin/env python3
"""Rebuild Part A Hill diversity/composition TIFF from saved analytical CSVs.

No bootstrap, permutation, selection or other analysis is rerun. The original
three-panel R plotting source is retained with display-only modifications.
"""
from pathlib import Path
import os
import re
from part_a_figure_presentation import patch_chunk

root = Path(__file__).resolve().parents[1]
snap = Path(os.environ.get("SEN_PARTA_SNAPSHOTS", str(root))).expanduser().resolve()
notebook = root / "workflow/A2_Genotype/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
out = root / "results/part_a/mdpi_label_corrections"
source_dir = snap / "results/part_a/full_master_original"
inputs = (
    source_dir / "QR_diversity_q0_q1_q2_bySource_with95CI.csv",
    source_dir / "QR_composition_bySource.csv",
)
for file in (notebook, *inputs):
    if not file.is_file() or file.stat().st_size == 0:
        raise SystemExit(f"[ERROR] Validated diversity input not found: {file}")
chunks = re.findall(r'(?ms)^\x60\x60\x60\{r[^\n}]*\}\s*\n(.*?)^\x60\x60\x60\s*$', notebook.read_text(encoding="utf-8"))
if len(chunks) != 30:
    raise SystemExit(f"[ERROR] Expected 30 notebook chunks, got {len(chunks)}")
original = chunks[21]
marker = "# G) Shannon + Simpson + composition patchwork"
if original.count(marker) != 1:
    raise SystemExit("[ERROR] Diversity plotting boundary missing")
start = original.rfind("# =========================================================", 0, original.index(marker))
if start < 0:
    raise SystemExit("[ERROR] Diversity section heading missing")
plotting = patch_chunk(22, original[start:])
header = r'''
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(readr)
})
root <- Sys.getenv("SEN_ROOT")
snap <- Sys.getenv("SEN_PARTA_SNAPSHOTS", unset = root)
out <- file.path(root, "results/part_a/mdpi_label_corrections")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "scripts/part_a_figure_label_helpers.R"))
source_dir <- file.path(snap, "results/part_a/full_master_original")
div_sum <- readr::read_csv(
  file.path(source_dir, "QR_diversity_q0_q1_q2_bySource_with95CI.csv"),
  show_col_types = FALSE)
comp_tbl <- readr::read_csv(
  file.path(source_dir, "QR_composition_bySource.csv"),
  show_col_types = FALSE)
source_levels <- c("Human", "Chicken-US", "Other-US", "International")
stopifnot(nrow(div_sum) == 4L,
          setequal(as.character(div_sum$Source), source_levels),
          setequal(as.character(comp_tbl$Source), source_levels),
          all(c("q1","q1_lo","q1_hi","q2","q2_lo","q2_hi") %in% names(div_sum)),
          all(c("marker_group","prevalence","p") %in% names(comp_tbl)))
marker_levels <- unique(as.character(comp_tbl$marker_group))
div_sum$Source <- factor(div_sum$Source, levels = source_levels)
comp_tbl$Source <- factor(comp_tbl$Source, levels = source_levels)
comp_tbl$marker_group <- factor(comp_tbl$marker_group, levels = marker_levels)
setwd(out)
'''
footer = '''
stopifnot(file.exists("QR_diversity_and_composition_side_by_side.tif"))
cat("[PASS] Three-panel diversity TIFF regenerated with source em dashes and italic gene legend.\\n")
'''
out.mkdir(parents=True, exist_ok=True)
target = out / "run_diversity_figure_only.R"
target.write_text(header + "\n" + plotting + "\n" + footer, encoding="utf-8")
print("[READY]", target)
