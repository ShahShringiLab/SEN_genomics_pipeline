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
echo "[1/4] Strict 33,831 x 51 parsed raw snapshot and metadata QC"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/00_parse_master_33831.R
echo "[AUDIT] Save raw malformed rows before/after parser and all corrected QR genotypes"
conda run --no-capture-output -n sen_stats \
  python scripts/audit_part_a2_csv_repairs.py
echo "[2/4] Prepare original Rmd with consistent temporal eligibility (all dated years <=2025)"
python3 scripts/prepare_original_part_a2.py "$NOTEBOOK"
echo "[3/4] Generate year × source × marker sampling/QC table and weights"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/07_selection_qc_weights_table.R
echo "[AUDIT] Validate temporal plot: 33,508 genomes; <=2000 bin has 38"
python3 scripts/audit_full_temporal_counts.py
echo "[4/4] Generate manuscript Table SS1: 2018-2025 by QR/FQR marker genotype, including QS/FQS"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/08_table_ss1_post2018_genotype.R
echo "[READY] analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
echo "[NOTE] Full-master Rmd statistics are unchanged. Only input paths, temporal eligibility/facet N, and Top-100 axis are adjusted."
echo "[NOTE] To execute the complete original notebook and QC followups, run bash scripts/run_part_a.sh."
