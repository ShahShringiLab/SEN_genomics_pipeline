#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
SRR_LIST="$REPO_ROOT/tests/smoke_test/smoke_srrs.txt"
KRAKEN_ENV="${SEN_SMOKE_KRAKEN_ENV:-sen_kraken}"
SEQSERO_ENV="${SEN_SMOKE_SEQSERO_ENV:-sen_seqsero2}"
SISTR_ENV="${SEN_SMOKE_SISTR_ENV:-sen_sistr}"
MLST_ENV="${SEN_SMOKE_MLST_ENV:-sen_mlst}"

KRAKEN_YAML="$REPO_ROOT/environments/02_kraken.yaml"
SEQSERO_YAML="$REPO_ROOT/environments/03_seqsero2.yaml"
SISTR_YAML="$REPO_ROOT/environments/03_sistr.yaml"
MLST_YAML="$REPO_ROOT/environments/04_mlst.yaml"

mkdir -p "$WORKDIR" "$WORKDIR/tmp" "$WORKDIR/tmp/kraken"
export SEN_ROOT="$WORKDIR"
export TMPDIR="$WORKDIR/tmp"
export SEN_KRAKEN_TMP="$WORKDIR/tmp/kraken"
export SEN_KRAKEN_CLEAN_DIR="$WORKDIR/Kraken_cleanup/clean_trimmed_fastq"

# shellcheck disable=SC1091
source "$REPO_ROOT/tests/smoke_test/checkpoint_lib.sh"

EXPECTED_SAMPLES="$(grep -cve '^[[:space:]]*$' "$SRR_LIST")"

echo "=================================================="
echo " SEN smoke test: Kraken -> SeqSero2/SISTR -> MLST"
echo "=================================================="
echo "[INFO] Repo:        $REPO_ROOT"
echo "[INFO] Workdir:     $WORKDIR"
echo "[INFO] Samples:     $EXPECTED_SAMPLES"
echo "[INFO] Checkpoints: $CHECKPOINT_ROOT"

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda is not available in PATH." >&2
  exit 1
}

create_env() {
  local env_name="$1" yaml="$2"
  echo "[INFO] Creating $env_name from $yaml ..."
  if conda env create --help 2>&1 | grep -q -- '--solver'; then
    conda env create --solver=libmamba -f "$yaml"
  else
    echo "[ERROR] This Conda installation does not expose --solver=libmamba." >&2
    exit 1
  fi
}

ensure_env() {
  local env_name="$1" yaml="$2" required_cmd="$3"
  if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
    create_env "$env_name" "$yaml"
  elif ! conda run -n "$env_name" command -v "$required_cmd" >/dev/null 2>&1; then
    echo "[WARN] $env_name exists but $required_cmd is unavailable; rebuilding."
    conda env remove -n "$env_name" -y
    create_env "$env_name" "$yaml"
  fi
}

validate_trim_inputs() {
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$WORKDIR/trimmed_fastq/${s}_trimmed_1.fastq.gz" ]] || return 1
    [[ -s "$WORKDIR/trimmed_fastq/${s}_trimmed_2.fastq.gz" ]] || return 1
  done < "$SRR_LIST"
}

KRAKEN_LOG="$WORKDIR/Kraken_cleanup/full_db_results.csv"
validate_kraken() {
  [[ -s "$KRAKEN_LOG" ]] || return 1
  [[ "$(( $(wc -l < "$KRAKEN_LOG") - 1 ))" -eq "$EXPECTED_SAMPLES" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$SEN_KRAKEN_CLEAN_DIR/${s}_pure_1.fastq.gz" ]] || return 1
    [[ -s "$SEN_KRAKEN_CLEAN_DIR/${s}_pure_2.fastq.gz" ]] || return 1
    [[ -s "$WORKDIR/Kraken_cleanup/kraken_reports/${s}.report" ]] || return 1
    awk -F',' -v id="$s" 'NR>1 && $1==id && $3>=30 {ok=1} END{exit !ok}' "$KRAKEN_LOG" || return 1
  done < "$SRR_LIST"
}

SEQSERO_SUMMARY="$WORKDIR/seqsero2_results/SeqSero2_summary.tsv"
validate_seqsero() {
  [[ -s "$SEQSERO_SUMMARY" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$WORKDIR/seqsero2_results/$s/SeqSero_result.tsv" ]] || return 1
  done < "$SRR_LIST"
  [[ "$(grep -ic 'Enteritidis' "$SEQSERO_SUMMARY" || true)" -ge "$EXPECTED_SAMPLES" ]]
}

SISTR_SUMMARY="$WORKDIR/sistr_results_run/sistr_master_summary.csv"
validate_sistr() {
  [[ -s "$SISTR_SUMMARY" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$WORKDIR/sistr_results_run/mini_assemblies/$s.fasta" ]] || return 1
    [[ -s "$WORKDIR/sistr_results_run/individual_csvs/$s.csv" ]] || return 1
  done < "$SRR_LIST"
  [[ "$(grep -ic 'Enteritidis' "$SISTR_SUMMARY" || true)" -ge "$EXPECTED_SAMPLES" ]]
}

MLST_SUMMARY="$WORKDIR/mlst_results_run/mlst_master_report.csv"
validate_mlst() {
  [[ -s "$MLST_SUMMARY" ]] || return 1
  [[ "$(( $(wc -l < "$MLST_SUMMARY") - 1 ))" -ge "$EXPECTED_SAMPLES" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    grep -q "$s" "$MLST_SUMMARY" || return 1
  done < "$SRR_LIST"
}

bootstrap_kraken_db() {
  if [[ -z "${SEN_KRAKEN_DB:-}" && -f "$REPO_ROOT/config/paths.env" ]]; then
    # shellcheck disable=SC1091
    source "$REPO_ROOT/config/paths.env"
  fi

  local default_db="$WORKDIR/databases/kraken2/k2_standard_20260626"
  if [[ -z "${SEN_KRAKEN_DB:-}" || \
        ! -s "${SEN_KRAKEN_DB:-}/hash.k2d" || \
        ! -s "${SEN_KRAKEN_DB:-}/opts.k2d" || \
        ! -s "${SEN_KRAKEN_DB:-}/taxo.k2d" ]]; then
    echo "[INFO] No complete Kraken2 database configured; bootstrapping pinned snapshot."
    conda run --no-capture-output -n "$KRAKEN_ENV" \
      env SEN_ROOT="$SEN_ROOT" TMPDIR="$TMPDIR" \
      bash "$REPO_ROOT/scripts/setup_kraken_db.sh"
    export SEN_KRAKEN_DB="$default_db"
  fi

  [[ -s "$SEN_KRAKEN_DB/hash.k2d" && -s "$SEN_KRAKEN_DB/opts.k2d" && -s "$SEN_KRAKEN_DB/taxo.k2d" ]] || {
    echo "[ERROR] Kraken2 database is incomplete after bootstrap." >&2
    exit 1
  }
  echo "[INFO] Kraken DB: $SEN_KRAKEN_DB"
}

validate_trim_inputs || {
  echo "[ERROR] Trimmed read inputs are incomplete. Run run_reads_qc.sh first." >&2
  exit 1
}

echo
echo "--------------------------------------------------"
echo "[STAGE] 1/4 Kraken2 filtering + Enterobacteriaceae extraction"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "kraken_filter" validate_kraken \
    "$REPO_ROOT/workflow/02_qc/06_kraken_cleanup.sh" "$KRAKEN_YAML" \
    "$REPO_ROOT/config/database_sources.env" "$SRR_LIST"; then
  ensure_env "$KRAKEN_ENV" "$KRAKEN_YAML" kraken2
  bootstrap_kraken_db
  conda run --no-capture-output -n "$KRAKEN_ENV" \
    env \
      SEN_ROOT="$SEN_ROOT" SEN_KRAKEN_DB="$SEN_KRAKEN_DB" \
      SEN_KRAKEN_TMP="$SEN_KRAKEN_TMP" \
      SEN_KRAKEN_PARALLEL_JOBS="${SEN_KRAKEN_PARALLEL_JOBS:-2}" \
      SEN_KRAKEN_THREADS="${SEN_KRAKEN_THREADS:-4}" \
      SEN_PIGZ_THREADS="${SEN_PIGZ_THREADS:-2}" \
      TMPDIR="$TMPDIR" \
      bash "$REPO_ROOT/workflow/02_qc/06_kraken_cleanup.sh"
  checkpoint_require_valid "kraken_filter" validate_kraken
  checkpoint_mark "kraken_filter" \
    "$REPO_ROOT/workflow/02_qc/06_kraken_cleanup.sh" "$KRAKEN_YAML" \
    "$REPO_ROOT/config/database_sources.env" "$SRR_LIST"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 2/4 SeqSero2"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "seqsero2" validate_seqsero \
    "$REPO_ROOT/workflow/03_typing/10_seqsero2.sh" "$SEQSERO_YAML" "$SRR_LIST"; then
  ensure_env "$SEQSERO_ENV" "$SEQSERO_YAML" SeqSero2_package.py
  conda run --no-capture-output -n "$SEQSERO_ENV" \
    env SEN_ROOT="$SEN_ROOT" SEN_SEQSERO2_JOBS="${SEN_SEQSERO2_JOBS:-2}" TMPDIR="$TMPDIR" \
    bash "$REPO_ROOT/workflow/03_typing/10_seqsero2.sh"
  checkpoint_require_valid "seqsero2" validate_seqsero
  checkpoint_mark "seqsero2" \
    "$REPO_ROOT/workflow/03_typing/10_seqsero2.sh" "$SEQSERO_YAML" "$SRR_LIST"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 3/4 SKESA + SISTR"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "sistr_skesa" validate_sistr \
    "$REPO_ROOT/workflow/03_typing/11_sistr.sh" "$SISTR_YAML" "$SRR_LIST"; then
  ensure_env "$SISTR_ENV" "$SISTR_YAML" skesa
  if ! conda run -n "$SISTR_ENV" command -v sistr >/dev/null 2>&1; then
    echo "[WARN] $SISTR_ENV is missing sistr; rebuilding."
    conda env remove -n "$SISTR_ENV" -y
    create_env "$SISTR_ENV" "$SISTR_YAML"
  fi
  conda run --no-capture-output -n "$SISTR_ENV" \
    env \
      SEN_ROOT="$SEN_ROOT" SEN_KRAKEN_CLEAN_DIR="$SEN_KRAKEN_CLEAN_DIR" \
      SEN_SISTR_JOBS="${SEN_SISTR_JOBS:-2}" \
      SEN_SKESA_CORES="${SEN_SKESA_CORES:-4}" \
      SEN_SKESA_MEMORY_GB="${SEN_SKESA_MEMORY_GB:-4}" \
      TMPDIR="$TMPDIR" \
      bash "$REPO_ROOT/workflow/03_typing/11_sistr.sh"
  checkpoint_require_valid "sistr_skesa" validate_sistr
  checkpoint_mark "sistr_skesa" \
    "$REPO_ROOT/workflow/03_typing/11_sistr.sh" "$SISTR_YAML" "$SRR_LIST"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 4/4 MLST"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "mlst" validate_mlst \
    "$REPO_ROOT/workflow/03_typing/12_mlst.sh" "$MLST_YAML" "$SRR_LIST"; then
  ensure_env "$MLST_ENV" "$MLST_YAML" mlst
  conda run --no-capture-output -n "$MLST_ENV" \
    env SEN_ROOT="$SEN_ROOT" SEN_MLST_JOBS="${SEN_MLST_JOBS:-4}" TMPDIR="$TMPDIR" \
    bash "$REPO_ROOT/workflow/03_typing/12_mlst.sh"
  checkpoint_require_valid "mlst" validate_mlst
  checkpoint_mark "mlst" \
    "$REPO_ROOT/workflow/03_typing/12_mlst.sh" "$MLST_YAML" "$SRR_LIST"
fi

post30="$(awk -F',' 'NR>1 && $3>=30 {n++} END{print n+0}' "$KRAKEN_LOG")"
pure_pairs="$(find "$SEN_KRAKEN_CLEAN_DIR" -maxdepth 1 -name '*_pure_1.fastq.gz' -size +0c | wc -l)"
seqsero_enteritidis="$(grep -ic 'Enteritidis' "$SEQSERO_SUMMARY" || true)"
sistr_enteritidis="$(grep -ic 'Enteritidis' "$SISTR_SUMMARY" || true)"
assembly_count="$(find "$WORKDIR/sistr_results_run/mini_assemblies" -maxdepth 1 -name '*.fasta' -size +0c | wc -l)"
mlst_rows="$(( $(wc -l < "$MLST_SUMMARY") - 1 ))"

echo
echo "=================================================="
echo " KRAKEN/TYPING SMOKE SUMMARY"
echo "=================================================="
cat "$KRAKEN_LOG"
echo
echo "Kraken-cleaned read pairs: $pure_pairs/$EXPECTED_SAMPLES"
echo "Post-Kraken depth >=30x:   $post30/$EXPECTED_SAMPLES"
echo "SeqSero2 Enteritidis hits: $seqsero_enteritidis"
echo "SISTR Enteritidis hits:    $sistr_enteritidis"
echo "SKESA assemblies:          $assembly_count/$EXPECTED_SAMPLES"
echo "MLST result rows:           $mlst_rows"
echo
echo "SMOKE KRAKEN/TYPING RESULT: PASS"
