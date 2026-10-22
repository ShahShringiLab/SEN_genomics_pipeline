#!/usr/bin/env python3
"""Reproduce the uploaded 30-chunk original Rmd on correct populations.

Original Rmd chunks 3-22: all 33,831 genomes.
Original Rmd chunks 24-30: the validated frozen 3,434 selected genomes.
Historical selection-generation chunk 23 is not rerun.
"""
from pathlib import Path
import hashlib
import re

root = Path(__file__).resolve().parents[1]
src_file = root / "analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
if not src_file.is_file():
    raise SystemExit("Run bash scripts/run_master_part_a2.sh first")
src = src_file.read_text(encoding="utf-8")
parts = re.findall(r'(?ms)^\x60\x60\x60\{r[^\n}]*\}\s*\n(.*?)^\x60\x60\x60\s*$', src)
assert len(parts) == 30, f"Expected complete original 30-chunk Rmd, got {len(parts)}"
assert "FINAL 3500-GENOME SELECTION SCRIPT" in parts[22]
assert "DOTPLOT BY YEAR" in parts[18]
assert "SEN_TableS1_final_source_decomposition.csv" in parts[16]

out = root / "results/part_a/master_33831/original_execution"
master_out = root / "results/part_a/full_master_original"
selected_out = root / "results/part_a/selected_original"
for d in (out, master_out, selected_out):
    d.mkdir(parents=True, exist_ok=True)

master = [
    'options(warn=1)',
    'suppressPackageStartupMessages({',
    'library(dplyr); library(tidyr); library(readr); library(stringr)',
    'library(ggplot2); library(lubridate); library(purrr); library(patchwork)',
    'library(scales); library(forcats)',
    '})',
    'setwd(file.path(Sys.getenv("SEN_ROOT"), "results/part_a/full_master_original"))',
    f'cat("[INFO] Original Rmd working-copy SHA256: {hashlib.sha256(src.encode()).hexdigest()}\\n")',
]
for index in range(3, 23):
    master.extend([f'\n# BEGIN USER ORIGINAL CHUNK {index}',
                   parts[index - 1],
                   f'# END USER ORIGINAL CHUNK {index}'])
    if index >= 4:
        master.extend([
            'if (!exists("SEN") || nrow(SEN) != 33831L || dplyr::n_distinct(SEN$Run) != 33831L)',
            f'  stop("33,831-genome population compromised after original chunk {index}")'])
    if index == 19:
        master.extend([
            'stopifnot(exists("den_year"), sum(den_year$N)==33508L)',
            'stopifnot(exists("src_N"), sum(src_N$N_total)==33508L)',
            'cat("[PASS] Temporal dot plot: 33,508 dated genomes; <=2000 included\\n")'])
master.append('cat("[PASS] All original full-master Rmd chunks 3-22: 33,831 genomes\\n")')
master_file = out / "run_original_A2_master.R"
master_file.write_text("\n".join(master)+"\n", encoding="utf-8")

canonical = 'file.path(Sys.getenv("SEN_ROOT"), "results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv")'
keys = 'file.path(Sys.getenv("SEN_ROOT"), "data/part_a/selection/Selected_updated_Isolate.csv")'
followups = [
    'setwd(file.path(Sys.getenv("SEN_ROOT"), "results/part_a/selected_original"))',
    'SEL <- readr::read_csv(' + canonical + ', show_col_types=FALSE)',
    'stopifnot(nrow(SEL)==3434L, dplyr::n_distinct(SEL$Isolate)==3434L)',
    'cat("[INFO] Selected analysis: 3,434 frozen QC-merged genomes, no resampling\\n")']
for index in range(24, 31):
    chunk = parts[index-1]
    chunk = chunk.replace('"Selected_3500_genomes.csv"', canonical)
    chunk = chunk.replace('"Finalupdated_SEN_Genomes1.csv"', canonical)
    chunk = chunk.replace('"Selected_updated_Isolate.csv"', keys)
    followups.extend([f'\n# BEGIN USER ORIGINAL CHUNK {index}',
                      chunk, f'# END USER ORIGINAL CHUNK {index}',
                      'if (!exists("SEN") || nrow(SEN)!=33831L) stop("Full master replaced by selected subset")',
                      'if (!exists("SEL") || nrow(SEL)!=3434L) stop("Selected cohort row count changed")'])
followups.append('cat("[PASS] Original selected followups (chunks 24-30): 3,434 QC-merged genomes\\n")')
all_file = out / "run_original_A2_all.R"
all_file.write_text(
    'source(file.path(Sys.getenv("SEN_ROOT"), "results/part_a/master_33831/original_execution/run_original_A2_master.R"))\n'
    + "\n".join(followups)+"\n", encoding="utf-8")
print(f"[READY] Original full-master (3-22): {master_file}")
print(f"[READY] Original full + selected followups (24-30): {all_file}")
print("[AUDIT] Original selection-generation chunk 23 intentionally skipped; existing QC-merged cohort reused")
