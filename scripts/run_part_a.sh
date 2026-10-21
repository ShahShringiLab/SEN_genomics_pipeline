#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SEN_ROOT="${SEN_ROOT:-$REPO_ROOT}"
cd "$SEN_ROOT"

echo "=================================================="
echo " SEN Part A: phenotype + genotype/cohort analysis"
echo "=================================================="

required=(
  data/part_a/phenotype/FDASEN_AMRPhenotype2.csv
  data/part_a/genotype/US-SEN_12192025.csv
  data/part_a/selection/Selected_updated_Isolate.csv
  data/part_a/selection/Finalupdated_SEN_Genomes1.csv
)

for f in "${required[@]}"; do
  [[ -s "$f" ]] || {
    echo "[FAIL] missing input: $f" >&2
    exit 1
  }
done

echo
echo "[VALIDATE] Part A inputs and Part A -> Part B handoff"
conda run --no-capture-output -n sen_stats   env SEN_ROOT="$SEN_ROOT"   python scripts/validate_part_a_inputs.py

run_r() {
  local script="$1"
  echo
  echo "[RUN] $script"
  conda run --no-capture-output -n sen_stats     env SEN_ROOT="$SEN_ROOT"     Rscript "$script"
}

echo
echo "[RUN] Master preparation, Rmd correction, and both year/QC sampling tables"
bash scripts/run_master_part_a2.sh
run_r analysis/part_a/01_phenotype_trends.R
run_r analysis/part_a/02_build_final_cohort.R
run_r analysis/part_a/03_genotype_marker_summary.R
run_r analysis/part_a/04_source_profiles.R
run_r analysis/part_a/05_diversity_composition.R

echo
echo "[VALIDATE] Generated Part A outputs"
conda run --no-capture-output -n sen_stats \
  env SEN_ROOT="$SEN_ROOT" \
  python scripts/validate_part_a_outputs.py

echo
echo "=================================================="
echo "[RUN] Original master Rmd chunks 4-21"
python3 scripts/build_original_a2_master_execution.py
MASTER_DIR="$SEN_ROOT/results/part_a/master_33831/original_execution"
if ! conda run --no-capture-output -n sen_stats env SEN_ROOT="$SEN_ROOT" \
  Rscript "$MASTER_DIR/run_original_A2_master.R" > "$MASTER_DIR/run.log" 2>&1; then
  echo "[FAIL] Original A2 execution; see $MASTER_DIR/run.log" >&2
  tail -n 60 "$MASTER_DIR/run.log" >&2
  exit 1
fi
grep -F "[PASS] Temporal dot plot:" "$MASTER_DIR/run.log" || {
  echo "[FAIL] Temporal dot plot audit did not appear in $MASTER_DIR/run.log" >&2
  exit 1
}
echo "[PASS] Original master analyses completed with validated 33,508-genome temporal denominator"
echo " SEN PART A: PASS"
echo "=================================================="
