#!/usr/bin/env bash
# Preserve original Part A2 order: all 33,831 first, selected cohort only afterward.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"
run_r() {
  conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" Rscript "$1"
}
echo "[1/3] Original 33,831 master: normalize, validate, export RAW 51-column XLSX/CSV."
run_r analysis/part_a/00_parse_master_33831.R
echo "[2/3] Original Part A2 master-population QR/FQR analyses and plots."
run_r analysis/part_a/06_master_33831_analysis.R
echo "[3/3] Downstream Part A1 and 3,434-genome selected-cohort A2; validate."
bash scripts/run_part_a.sh
echo "[PASS] Master 33,831 analyses preceded downstream 3,434 cohort analyses."
