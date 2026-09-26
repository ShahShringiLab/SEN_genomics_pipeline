#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

ENV_NAME="${SEN_SMOKE_READS_ENV:-sen_reads_qc}"
ENV_YAML="$REPO_ROOT/environments/01_reads_qc.yaml"
WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
SRR_LIST="$REPO_ROOT/tests/smoke_test/smoke_srrs.txt"

mkdir -p "$WORKDIR" "$WORKDIR/tmp"
export SEN_ROOT="$WORKDIR"
export SEN_SRR_LIST="$SRR_LIST"
export TMPDIR="$WORKDIR/tmp"
export SEN_TMP_DIR="$WORKDIR/tmp"

# shellcheck disable=SC1091
source "$REPO_ROOT/tests/smoke_test/checkpoint_lib.sh"

EXPECTED_SAMPLES="$(grep -cve '^[[:space:]]*$' "$SRR_LIST")"
EXPECTED_FASTQ=$((EXPECTED_SAMPLES * 2))

echo "=================================================="
echo " SEN smoke test: download -> fastp -> QC -> coverage"
echo "=================================================="
echo "[INFO] Repo:        $REPO_ROOT"
echo "[INFO] Workdir:     $WORKDIR"
echo "[INFO] Samples:     $EXPECTED_SAMPLES"
echo "[INFO] Checkpoints: $CHECKPOINT_ROOT"
echo "[INFO] TMPDIR:      $TMPDIR"

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda is not available in PATH." >&2
  exit 1
}

create_reads_env() {
  echo "[INFO] Creating $ENV_NAME from $ENV_YAML ..."
  if conda env create --help 2>&1 | grep -q -- '--solver'; then
    conda env create --solver=libmamba -f "$ENV_YAML"
  else
    echo "[ERROR] This Conda installation does not expose --solver=libmamba." >&2
    exit 1
  fi
}

ensure_reads_env() {
  if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    create_reads_env
  elif ! conda run -n "$ENV_NAME" python -c 'import pkg_resources' >/dev/null 2>&1 || \
       ! conda run -n "$ENV_NAME" multiqc --version >/dev/null 2>&1; then
    echo "[WARN] Existing $ENV_NAME is incompatible; rebuilding."
    conda env remove -n "$ENV_NAME" -y
    create_reads_env
  fi
}

sample_file_exists() {
  local root="$1" pattern="$2"
  find "$root" -type f -name "$pattern" -size +0c -print -quit 2>/dev/null | grep -q .
}

validate_download() {
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    sample_file_exists "$WORKDIR/fastq" "${s}_1.fastq.gz" || return 1
    sample_file_exists "$WORKDIR/fastq" "${s}_2.fastq.gz" || return 1
  done < "$SRR_LIST"
}

validate_raw_qc() {
  local n
  n="$(find "$WORKDIR/raw_fastq_qc" -maxdepth 1 -type f -name '*_fastqc.zip' -size +0c 2>/dev/null | wc -l)"
  [[ "$n" -eq "$EXPECTED_FASTQ" && -s "$WORKDIR/raw_fastq_qc/multiqc_report.html" ]]
}

validate_trim() {
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$WORKDIR/trimmed_fastq/${s}_trimmed_1.fastq.gz" ]] || return 1
    [[ -s "$WORKDIR/trimmed_fastq/${s}_trimmed_2.fastq.gz" ]] || return 1
    [[ -s "$WORKDIR/logs/${s}.json" ]] || return 1
  done < "$SRR_LIST"
}

validate_trimmed_qc() {
  local n
  n="$(find "$WORKDIR/trimmed_fastq_qc" -maxdepth 1 -type f -name '*_fastqc.zip' -size +0c 2>/dev/null | wc -l)"
  [[ "$n" -eq "$EXPECTED_FASTQ" && -s "$WORKDIR/trimmed_fastq_qc/trimmed_data_multiqc_report.html" ]]
}

validate_coverage() {
  local csv="$WORKDIR/trimmed_fastq_qc/all_samples_coverage.csv"
  [[ -s "$csv" ]] || return 1
  [[ "$(( $(wc -l < "$csv") - 1 ))" -eq "$EXPECTED_SAMPLES" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    awk -F',' -v id="$s" 'NR>1 && $1==id {found=1} END{exit !found}' "$csv" || return 1
  done < "$SRR_LIST"
}

run_in_reads_env() {
  local script="$1"
  conda run --no-capture-output -n "$ENV_NAME" \
    env \
      SEN_ROOT="$SEN_ROOT" \
      SEN_SRR_LIST="$SEN_SRR_LIST" \
      SEN_TMP_DIR="$SEN_TMP_DIR" \
      SEN_FASTP_JOBS="${SEN_FASTP_JOBS:-2}" \
      SEN_FASTP_THREADS_PER_JOB="${SEN_FASTP_THREADS_PER_JOB:-4}" \
      SEN_FASTQC_RAW_JOBS="${SEN_FASTQC_RAW_JOBS:-2}" \
      SEN_FASTQC_TRIMMED_JOBS="${SEN_FASTQC_TRIMMED_JOBS:-2}" \
      SEN_FASTQC_THREADS_PER_JOB="${SEN_FASTQC_THREADS_PER_JOB:-2}" \
      SEN_COVERAGE_JOBS="${SEN_COVERAGE_JOBS:-2}" \
      TMPDIR="$TMPDIR" \
      bash "$script"
}

run_stage() {
  local number="$1" name="$2" validator="$3" script="$4"
  shift 4
  local sigfiles=("$@")

  echo
  echo "--------------------------------------------------"
  echo "[STAGE] $number $name"
  echo "--------------------------------------------------"

  if checkpoint_should_skip "$name" "$validator" "${sigfiles[@]}"; then
    return 0
  fi

  ensure_reads_env
  run_in_reads_env "$script"
  checkpoint_require_valid "$name" "$validator"
  checkpoint_mark "$name" "${sigfiles[@]}"
}

run_stage "1/5" "reads_download" validate_download \
  "$REPO_ROOT/workflow/01_reads/01_download_fastq.sh" \
  "$REPO_ROOT/workflow/01_reads/01_download_fastq.sh" "$ENV_YAML" "$SRR_LIST"

run_stage "2/5" "raw_fastqc" validate_raw_qc \
  "$REPO_ROOT/workflow/02_qc/02_fastqc_raw.sh" \
  "$REPO_ROOT/workflow/02_qc/02_fastqc_raw.sh" "$ENV_YAML" "$SRR_LIST"

run_stage "3/5" "fastp_trim" validate_trim \
  "$REPO_ROOT/workflow/03_typing/09_trim_fastp.sh" \
  "$REPO_ROOT/workflow/03_typing/09_trim_fastp.sh" "$ENV_YAML" "$SRR_LIST"

run_stage "4/5" "trimmed_fastqc" validate_trimmed_qc \
  "$REPO_ROOT/workflow/02_qc/03_fastqc_trimmed.sh" \
  "$REPO_ROOT/workflow/02_qc/03_fastqc_trimmed.sh" "$ENV_YAML" "$SRR_LIST"

run_stage "5/5" "coverage_estimate" validate_coverage \
  "$REPO_ROOT/workflow/02_qc/04_estimate_coverage.sh" \
  "$REPO_ROOT/workflow/02_qc/04_estimate_coverage.sh" "$ENV_YAML" "$SRR_LIST"

raw_count="$(find "$WORKDIR/fastq" -type f -name '*.fastq.gz' -size +0c 2>/dev/null | wc -l)"
trim_count="$(find "$WORKDIR/trimmed_fastq" -type f -name '*_trimmed_*.fastq.gz' -size +0c 2>/dev/null | wc -l)"

echo
echo "=================================================="
echo " SMOKE READ/QC SUMMARY"
echo "=================================================="
echo "Expected paired FASTQs: $EXPECTED_FASTQ"
echo "Raw FASTQs found:       $raw_count"
echo "Trimmed FASTQs found:   $trim_count"
echo "Coverage table:         $WORKDIR/trimmed_fastq_qc/all_samples_coverage.csv"
echo
cat "$WORKDIR/trimmed_fastq_qc/all_samples_coverage.csv"
echo
echo "SMOKE READ/QC RESULT: PASS"
