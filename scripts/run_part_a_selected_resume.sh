#!/usr/bin/env bash
# Resume at original selected-genome analyses after full 33,831-genome Rmd completed.
# Reuses the existing master raw+Rmd-derived snapshot; no master statistics rerun.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SEN_ROOT="$ROOT"
cd "$ROOT"

RAW="results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv"
DERIVED="results/part_a/full_master_original/SEN_33831_original_Rmd_derived.csv"
QC="results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv"
for f in "$RAW" "$DERIVED" "$QC"; do
  if [[ ! -s "$f" ]]; then
    echo "[FAIL] Saved prerequisite missing: $f" >&2
    echo "[ACTION] Run bash scripts/run_part_a.sh for a full original rerun." >&2
    exit 1
  fi
done

echo "[1/5] Recheck original Rmd integrity and rebuild updated execution scripts"
python3 scripts/check_original_rmd_contract.py

WORK="results/part_a/master_33831/original_execution"
echo "[2/5] Resume original selected chunks 24-30, without recomputing full master"
if ! conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript "$WORK/run_original_A2_selected_resume.R" \
  > "$WORK/selected_resume.log" 2>&1; then
  echo "[FAIL] Original selected analysis. See $WORK/selected_resume.log" >&2
  tail -n 90 "$WORK/selected_resume.log" >&2
  exit 1
fi

echo "[3/5] Validate original full-master results and selected QC"
python3 scripts/validate_full_master_repro.py

echo "[ADD-ON] International decomposed animal-source chart (N=71)"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/10_international_animal_source_plot.R

echo "[4/5] Generate matching Top-5 charts and enriched 33,831/3,434 metadata"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/09_temporal_top5_and_enriched_metadata.R

echo "[5/5] Validate enriched outputs and plot denominator contracts"
python3 scripts/validate_enriched_metadata_outputs.py
grep -F "[PASS] Original selected followups" "$WORK/selected_resume.log"
echo "[PASS] Selected original chunks plus new Top-5/master+QC metadata completed."
echo "[RESULTS] results/part_a/reconciled_metadata/"
echo "[RESULTS] results/part_a/temporal_top5/"
