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
echo " SEN PART A: PASS"
echo "=================================================="
