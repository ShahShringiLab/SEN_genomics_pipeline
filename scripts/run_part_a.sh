#!/usr/bin/env bash
# Original SEN Part A: separate full master from validated selected cohort.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SEN_ROOT="$ROOT"
cd "$ROOT"

echo "=============================================================="
echo " Original SEN Part A: 33,831 full master; 3,434 selected for QC"
echo "=============================================================="
required=(
  data/part_a/phenotype/FDASEN_AMRPhenotype2.csv
  data/part_a/genotype/US-SEN_12192025.csv
  data/part_a/selection/Selected_updated_Isolate.csv
  data/part_a/selection/Finalupdated_SEN_Genomes1.csv
  data/part_a/selection/Selected_US-SEN_12192025_QCAUDIT_R1_qr_clade.csv
  provenance/SEN_Genomes2_uploaded_original_20261021.Rmd
)
for file in "${required[@]}"; do
  [[ -s "$file" ]] || { echo "[FAIL] Missing original input: $file" >&2; exit 1; }
done

echo "[PREFLIGHT] Check exact original 30-chunk Rmd and dataset boundaries"
python3 scripts/check_original_rmd_contract.py

echo "[1/6] Validate original datasets and frozen selection keys"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  python scripts/validate_part_a_inputs.py

echo "[2/6] Correct master CSV loading; prepare original 30-chunk Rmd; sampling audits"
bash scripts/run_master_part_a2.sh

echo "[3/6] Reconcile fixed 3,434 selected metadata with core genotypes and clades"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/02_build_final_cohort.R

echo "[4/6] Merge verified QC edits onto complete 3,434 metadata (never 33,831)"
python3 scripts/merge_qc_full_metadata.py

echo "[5/6] Extract original analyses and separate output populations"
python3 scripts/build_original_a2_master_execution.py

MASTER_DIR="$ROOT/results/part_a/master_33831/original_execution"
echo "[6/6] Execute original analyses: 33,831 full master; 3,434 selected followups"
if ! conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
    Rscript "$MASTER_DIR/run_original_A2_all.R" > "$MASTER_DIR/run.log" 2>&1; then
  echo "[FAIL] Original Rmd execution. Review $MASTER_DIR/run.log" >&2
  tail -n 90 "$MASTER_DIR/run.log" >&2
  exit 1
fi

python3 scripts/validate_full_master_repro.py

echo "[ADD-ON] Seven-bucket Top-5 temporal charts and reconciled master/selected metadata"
conda run --no-capture-output -n sen_stats env SEN_ROOT="$ROOT" \
  Rscript analysis/part_a/09_temporal_top5_and_enriched_metadata.R
python3 scripts/validate_enriched_metadata_outputs.py

grep -F "[PASS] Temporal dot plot:" "$MASTER_DIR/run.log"
grep -F "[PASS] All original full-master Rmd chunks" "$MASTER_DIR/run.log"
grep -F "[PASS] Original selected followups" "$MASTER_DIR/run.log"

echo "[PASS] Original Part A analysis and QC completed."
echo "[RESULTS] 33,831 master figures: results/part_a/full_master_original/"
echo "[RESULTS] Original 3,434 analyses: results/part_a/selected_original/"
echo "[RESULTS] 3,434 metadata/QC: results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv"
echo "[RESULTS] Master/selected enriched metadata: results/part_a/reconciled_metadata/"
echo "[RESULTS] Seven-bucket full/selected Top-5 temporal stacks: results/part_a/temporal_top5/"
