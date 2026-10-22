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
    chunk = parts[index - 1]
    if index == 18:
        # Original chunk 18 used obsolete Source_group labels
        # "International-Non-US" and "Other-US-US"; chunk 17 already
        # builds the validated Final_source values International/Other-US.
        # Correct ONLY this impossible input filter/recode in the generated
        # execution script. Preserve the animal and QR grouping algorithm.
        obsolete = '''    Source_group %in% c(
      "International-Non-US",
      "Other-US-US"
    )'''
        actual = '''    Final_source %in% c(
      "International",
      "Other-US"
    )'''
        obsolete_recode = '''    Group = recode(
      Source_group,
      "International-Non-US" = "International",
      "Other-US-US" = "Other-US"
    ),'''
        actual_recode = '''    Group = Final_source,'''
        assert chunk.count(obsolete) == 1 and chunk.count(obsolete_recode) == 1
        chunk = chunk.replace(obsolete, actual, 1).replace(obsolete_recode, actual_recode, 1)
    master.extend([f'\n# BEGIN USER ORIGINAL CHUNK {index}',
                   chunk,
                   f'# END USER ORIGINAL CHUNK {index}'])
    if index >= 4:
        master.extend([
            'if (!exists("SEN") || nrow(SEN) != 33831L || dplyr::n_distinct(SEN$Run) != 33831L)',
            f'  stop("33,831-genome population compromised after original chunk {index}")'])
    if index == 17:
        # Reporting-only audit: preserve exact original source assignment and
        # allow traceability of the 22-genome Swine/Others manuscript difference.
        master.extend([
            'source_audit <- SEN %>% dplyr::filter(Final_source %in% c("Other-US", "International")) %>%',
            '  dplyr::transmute(Run, Isolate, \x60Isolation source\x60, Source_group, Final_source,',
            '    Primary_source, Other_source_group, Isolate_identifiers=\x60Isolate identifiers\x60,',
            '    Food_origin=\x60Food origin\x60)',
            'stopifnot(nrow(source_audit)==653L)',
            'readr::write_csv(source_audit, "SEN_fullmaster_source_isolate_audit.csv")',
            'swine_audit <- source_audit %>% dplyr::filter(Final_source=="Other-US", Other_source_group=="Swine") %>%',
            '  dplyr::mutate(matched_term=stringr::str_extract(stringr::str_to_lower(\x60Isolation source\x60),',
            '    "\\\\b(swine|sow|pork|porcine)\\\\b"))',
            'readr::write_csv(swine_audit, "SEN_OtherUS_swine_isolate_audit.csv")',
            'cat("[AUDIT] Full-master source-isolate and Swine-category audit exported\\n")',
        ])
    if index == 22:
        # Export the original Rmd-derived metadata before any selected analysis
        # has an opportunity to mutate SEN. Keep all 51 raw fields separately.
        master.extend([
            'raw_header <- readr::read_csv(file.path(Sys.getenv("SEN_ROOT"),',
            '  "results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv"),',
            '  n_max=0, show_col_types=FALSE)',
            'new_derived <- setdiff(names(SEN), names(raw_header))',
            'stopifnot(length(new_derived)>0, all(c("Run","Isolate","Final_source",',
            '  "Primary_source","Other_source_group") %in% names(SEN)))',
            'rmd_export <- SEN[, c("Run","Isolate",new_derived), drop=FALSE]',
            'rmd_export <- rmd_export %>% dplyr::mutate(dplyr::across(',
            '  where(is.list), ~vapply(.x,function(x) paste(as.character(unlist(x)),collapse=";"),character(1))))',
            'stopifnot(nrow(rmd_export)==33831L,dplyr::n_distinct(rmd_export$Run)==33831L)',
            'readr::write_csv(rmd_export,"SEN_33831_original_Rmd_derived.csv",na="")',
            'cat("[PASS] Original Rmd-derived full-master metadata exported: ",',
            '  length(new_derived)," fields for 33,831 genomes\\n",sep="")',
        ])
    if index == 19:
        master.extend([
            'stopifnot(exists("den_year"), sum(den_year$N)==33508L)',
            'stopifnot(exists("src_N"), all(c("Source", "Source_lab") %in% names(src_N)))',
            'stopifnot(dplyr::n_distinct(src_N$Source)==4L)',
            'stopifnot(sum(!is.na(SEN$Collection_year) & as.integer(SEN$Collection_year)<=2025L)==33508L)',
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
    'cat("[INFO] Selected analysis: 3,434 frozen QC-merged genomes, no resampling\\n")',
    '# Original chunk 25 expects SEN_plot/SEL_plot created in earlier interactive',
    '# sessions. Reconstruct ONLY these input views, using original Period rules.',
    '# Keep master SEN at 33,831 rows; do not resample or change original tests.',
    'stopifnot("FQR_genotype" %in% names(SEN))',
    '# Original stochastic chunk 23 (not rerun) also created the',
    '# QR_genotype alias for later Rmd chunks 25 and 30. Rehydrate',
    '# its binary marker alias without touching master Rmd results.',
    'SEN$QR_genotype <- as.integer(SEN$FQR_genotype)',
    'stopifnot(length(SEN$QR_genotype)==33831L,',
    '  all(SEN$QR_genotype %in% c(0L,1L)))',
    'make_period_original <- function(year) dplyr::case_when(',
    '  year < 2018L ~ "Pre-2018",',
    '  year >= 2018L & year <= 2020L ~ "2018-2020",',
    '  year >= 2021L & year <= 2025L ~ "2021-2025",',
    '  TRUE ~ NA_character_)',
    'SEN_plot <- SEN %>% dplyr::mutate(row_id=dplyr::row_number(),',
    '  Collection_year=as.integer(Collection_year),',
    '  QR_genotype=as.integer(FQR_genotype),',
    '  QR_status=dplyr::if_else(QR_genotype==1L,"QR/FQR marker","QS/FQS"),',
    '  Period=make_period_original(Collection_year)) %>%',
    '  dplyr::filter(!is.na(Period))',
    'SEL_plot <- SEL %>% dplyr::mutate(Collection_year=as.integer(Collection_year),',
    '  QR_genotype=as.integer(QR_genotype),',
    '  QR_status=dplyr::if_else(QR_genotype==1L,"QR/FQR marker","QS/FQS"),',
    '  Period=make_period_original(Collection_year)) %>%',
    '  dplyr::filter(!is.na(Period))',
    'stopifnot(nrow(SEN)==33831L, nrow(SEN_plot)==33508L,',
    '  nrow(SEL)==3434L, nrow(SEL_plot)==3434L,',
    '  all(c("Period","QR_genotype","QR_status","row_id") %in% names(SEN_plot)),',
    '  all(c("Period","QR_genotype","QR_status","row_id") %in% names(SEL_plot)),',
    '  sum(SEN_plot$QR_genotype==1L)+sum(SEN_plot$QR_genotype==0L)==33508L,',
    '  sum(SEL_plot$QR_genotype==1L)==2288L,',
    '  sum(SEL_plot$QR_genotype==0L)==1146L)',
    'cat("[PASS] Original interactive period views restored: SEN_plot 33,508 dated, SEL_plot 3,434 selected\\n")']
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
