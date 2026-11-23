#!/usr/bin/env python3
"""Build a figure-only R runner from the original validated Part A2 source.

This runs original Rmd chunk 6 (Figure 2 only) and chunks 19, 20 and 21
(Figure 3 and the four separately saved S2/S3 panels), with display-only
patches. Original Rmd, table outputs, and prior manuscript figures are intact.
"""
from __future__ import annotations

import re
from pathlib import Path
from part_a_figure_presentation import patch_chunk

root = Path(__file__).resolve().parents[1]
prepared = root / "workflow/A2_Genotype/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
snapshot_root = Path(__import__("os").environ.get("SEN_PARTA_SNAPSHOTS", str(root))).expanduser().resolve()
raw = snapshot_root / "results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv"
derived = snapshot_root / "results/part_a/full_master_original/SEN_33831_original_Rmd_derived.csv"
out = root / "results/part_a/mdpi_label_corrections"
out.mkdir(parents=True, exist_ok=True)

for path in (prepared, raw, derived):
    if not path.is_file():
        raise SystemExit(
            f"[MISSING] {path}\n"
            "First complete the existing Part A workflow with "
            "'bash scripts/run_part_a.sh'; no analytical inputs will be replaced."
        )

source = prepared.read_text(encoding="utf-8")
chunks = re.findall(r'(?ms)^\x60\x60\x60\{r[^\n}]*\}\s*\n(.*?)^\x60\x60\x60\s*$', source)
if len(chunks) != 30:
    raise SystemExit(f"[ERROR] Expected original 30 code chunks, found {len(chunks)}")

# The original figure 2 chunk generates additional analytical tables and
# unrelated plots after its first three file exports; do not execute those.
fig2 = chunks[5]
boundary = "# =========================\n# 5) Per-isolate marker counts"
if fig2.count(boundary) != 1:
    raise SystemExit("[ERROR] Figure 2-only execution boundary not uniquely found")
fig2 = fig2.split(boundary, 1)[0]

header = r'''
options(warn = 1)
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(stringr)
  library(ggplot2); library(lubridate); library(purrr); library(scales)
  library(forcats)
})
root <- Sys.getenv("SEN_ROOT")
if (!nzchar(root)) stop("SEN_ROOT must be provided.")
snapshots_root <- Sys.getenv("SEN_PARTA_SNAPSHOTS", unset = root)
out <- file.path(root, "results/part_a/mdpi_label_corrections")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "scripts/part_a_figure_label_helpers.R"))
ggsave <- sen_ggsave_part_a
# Fail before generating files if the requested typography cannot be parsed.
stopifnot(identical(unname(sen_pctbin_labels(c("1-5", "20-40"))),
                    c("1\u20145", "20\u201440")))
stopifnot(identical(unname(sen_marker_plotmath(c("Q/FQ-susceptible"))),
                    c('"QS/FQS"')))
stopifnot(identical(unname(sen_source_plotmath("Human (N=26296)")),
                    'paste("Human (", italic(n), " = 26,296)")'))
invisible(parse(text = sen_marker_plotmath(c("gyrA_D87Y", "qnrB19"))))
invisible(parse(text = sen_source_plotmath("Chicken-US (N=6588)")))

master_raw <- readr::read_csv(
  file.path(snapshots_root, "results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)
master_derived <- readr::read_csv(
  file.path(snapshots_root, "results/part_a/full_master_original/SEN_33831_original_Rmd_derived.csv"),
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE
)
stopifnot(nrow(master_raw) == 33831L, nrow(master_derived) == 33831L,
          !anyDuplicated(master_raw$Run), !anyDuplicated(master_derived$Run),
          all(c("Run", "Isolate") %in% names(master_derived)),
          !any(names(master_derived)[!names(master_derived) %in%
               c("Run", "Isolate")] %in% names(master_raw)))
SEN <- dplyr::left_join(master_raw, master_derived, by = c("Run", "Isolate"))
stopifnot(nrow(SEN) == 33831L, dplyr::n_distinct(SEN$Run) == 33831L)
setwd(out)
cat("[INFO] Regenerating manuscript Figures 2, 3, S2 and S3 from frozen 33,831-genome metadata.\n")
'''

code = [header]
for index, original_chunk in (
    (6, fig2), (19, chunks[18]), (20, chunks[19]), (21, chunks[20])
):
    chunk = patch_chunk(index, original_chunk)
    code += [
        f'\n# ===== ORIGINAL FULL-MASTER CHUNK {index} (display-only patch) =====\n',
        chunk,
        f'\ncat("[PASS] Figure source chunk {index} complete\\n")\n',
    ]
code.append('''
stopifnot(
  all(file.exists(file.path(out, c(
    "FQR_marker_counts_top100.tif",
    "FQ_burden.tif",
    "FQ_burden_5yr.tif",
    "FQ_QRFQR_vs_QFQsusceptible_PERCENT_by_year.tif",
    "FQ_QRFQR_vs_QFQsusceptible_COUNT_by_year.tif",
    "FQ_QSFQS_TOP5_QRFQR_PERCENT_by_year.tif",
    "FQ_QSFQS_TOP5_QRFQR_COUNT_by_year.tif"
  ))))
)
stopifnot(all(c("Human (N=26296)", "Chicken-US (N=6588)",
                  "Other-US (N=553)", "International (N=71)") %in%
              as.character(src_N$Source_lab)))
cat("[PASS] Figure 2, Figure 3, the 5-year companion, and four independent S2/S3 files regenerated.\\n")
''')
script = out / "run_figures_2_3_s2_s3.R"
script.write_text("\n".join(code) + "\n", encoding="utf-8")
print(f"[READY] Standalone publication figure execution script: {script}")
print("[AUDIT] Source Rmd and all statistical computations retained unchanged.")
print("[OUTPUT] " + str(out))
