#!/usr/bin/env bash
# Master Part A2: source CSV -> corrected 51-column XLSX -> original Part A -> 33K analysis.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"
echo "[1/3] Normalize and audit 33,831-record master; export UNMODIFIED 51-column XLSX."
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
 Rscript analysis/part_a/00_parse_master_33831.R
echo "[2/3] Run and validate existing Part A1 and A2 selected-cohort analyses."
bash scripts/run_part_a.sh
echo "[3/3] Produce original-population Part A2 derivatives, audits, and plots."
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
 Rscript analysis/part_a/06_master_33831_analysis.R
echo "[PASS] Raw 51-column XLSX + original Part A + master Part A2 complete."
