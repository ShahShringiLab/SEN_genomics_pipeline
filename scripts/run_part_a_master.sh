#!/usr/bin/env bash
# Run the corrected 33,831-record raw export FIRST, then the validated
# phenotype/selected analyses, then original-population Part A2 analyses.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"

if ! command -v conda >/dev/null 2>&1; then
  if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
  elif [[ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
  else
    echo "[ERROR] Activate Conda before running." >&2
    exit 1
  fi
fi

if ! conda run -n sen_stats Rscript -e 'stopifnot(requireNamespace("writexl", quietly=TRUE))'; then
  echo "[ERROR] Install writexl into sen_stats:" >&2
  echo "  conda install -n sen_stats -c conda-forge r-writexl -y" >&2
  exit 1
fi

echo "[RUN] Stage 0: 33,831-row corrected parse; 51 raw columns; XLSX + CSV"
conda run --no-capture-output -n sen_stats \
  env SEN_ROOT="$ROOT" Rscript analysis/part_a/00_parse_master_33831.R

echo "[RUN] Existing Part A1 and selected-cohort A2 + input/output validation"
bash scripts/run_part_a.sh

echo "[RUN] Master-population original Part A2: 33,831 rows"
conda run --no-capture-output -n sen_stats \
  env SEN_ROOT="$ROOT" Rscript analysis/part_a/06_master_33831_analysis.R

echo "[PASS] MASTER PART A2 + EXISTING PART A VALIDATED"
echo "Raw XLSX: results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.xlsx"
echo "Derived XLSX: results/part_a/master_33831/SEN_33831_DERIVED_ANALYSIS_READY.xlsx"
echo "Results: results/part_a/master_33831/"
