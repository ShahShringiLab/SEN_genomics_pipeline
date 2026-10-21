#!/usr/bin/env bash
# Original 33,831-genome Part A2, preserve original uploaded Rmd sequence/figures.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"
NOTEBOOK="${1:-$ROOT/provenance/SEN_Genomes2_uploaded_original_20261021.Rmd}"
if [[ ! -f "$NOTEBOOK" ]]; then
  echo "[ERROR] Exact uploaded original notebook not found: $NOTEBOOK" >&2
  echo "Copy your provided SEN_Genomes2.Rmd there or pass its path as the first argument." >&2
  exit 2
fi
echo "[1/3] Strict 33,831 x 51 parsed raw snapshot and metadata QC"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/00_parse_master_33831.R
echo "[2/3] Prepare original Rmd with parsed master and corrected temporal facet N"
python3 scripts/prepare_original_part_a2.py "$NOTEBOOK"
echo "[3/3] Generate year × source × marker sampling/QC table and weights"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/07_selection_qc_weights_table.R
echo "[READY] analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
echo "[NOTE] This preserves exact notebook analyses and plots; it does NOT silently run alternative analyses."
echo "[NOTE] Render/run the notebook only after resolving its original local input/output paths and exploratory state."
