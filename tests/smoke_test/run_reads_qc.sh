#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

ENV_NAME="${SEN_SMOKE_READS_ENV:-sen_reads_qc}"
ENV_YAML="$REPO_ROOT/environments/01_reads_qc.yaml"
WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
SRR_LIST="$REPO_ROOT/tests/smoke_test/smoke_srrs.txt"

echo "=================================================="
echo " SEN smoke test: download -> fastp -> QC -> coverage"
echo "=================================================="
echo "[INFO] Repo:     $REPO_ROOT"
echo "[INFO] Workdir:  $WORKDIR"
echo "[INFO] SRR list: $SRR_LIST"
echo "[INFO] Conda env:$ENV_NAME"

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda is not available in PATH." >&2
  exit 1
}

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[INFO] Creating $ENV_NAME from $ENV_YAML ..."

  echo "[INFO] Using conda libmamba solver."
  if conda env create --help 2>&1 | grep -q -- '--solver'; then
    conda env create --solver=libmamba -f "$ENV_YAML"
  else
    echo "[ERROR] This Conda installation does not expose --solver=libmamba." >&2
    echo "[ERROR] Please update Conda or install conda-libmamba-solver in base." >&2
    exit 1
  fi
fi

mkdir -p "$WORKDIR" "$WORKDIR/tmp"

export SEN_ROOT="$WORKDIR"
export SEN_SRR_LIST="$SRR_LIST"
export TMPDIR="$WORKDIR/tmp"
export SEN_TMP_DIR="$WORKDIR/tmp"

echo "[INFO] TMPDIR:   $TMPDIR"

# Laptop-safe defaults; override before launch if desired.
export SEN_PREFETCH_JOBS="${SEN_PREFETCH_JOBS:-2}"
export SEN_DUMP_JOBS="${SEN_DUMP_JOBS:-2}"
export SEN_THREADS_PER_DUMP="${SEN_THREADS_PER_DUMP:-4}"
export SEN_FASTP_JOBS="${SEN_FASTP_JOBS:-2}"
export SEN_FASTP_THREADS_PER_JOB="${SEN_FASTP_THREADS_PER_JOB:-4}"
export SEN_FASTQC_RAW_JOBS="${SEN_FASTQC_RAW_JOBS:-2}"
export SEN_FASTQC_TRIMMED_JOBS="${SEN_FASTQC_TRIMMED_JOBS:-2}"
export SEN_FASTQC_THREADS_PER_JOB="${SEN_FASTQC_THREADS_PER_JOB:-2}"
export SEN_COVERAGE_JOBS="${SEN_COVERAGE_JOBS:-2}"

run_stage() {
  local label="$1"
  local script="$2"
  echo
  echo "--------------------------------------------------"
  echo "[STAGE] $label"
  echo "--------------------------------------------------"
  conda run --no-capture-output -n "$ENV_NAME"     env       SEN_ROOT="$SEN_ROOT"       SEN_SRR_LIST="$SEN_SRR_LIST"       SEN_PREFETCH_JOBS="$SEN_PREFETCH_JOBS"       SEN_DUMP_JOBS="$SEN_DUMP_JOBS"       SEN_THREADS_PER_DUMP="$SEN_THREADS_PER_DUMP"       SEN_FASTP_JOBS="$SEN_FASTP_JOBS"       SEN_FASTP_THREADS_PER_JOB="$SEN_FASTP_THREADS_PER_JOB"       SEN_FASTQC_RAW_JOBS="$SEN_FASTQC_RAW_JOBS"       SEN_FASTQC_TRIMMED_JOBS="$SEN_FASTQC_TRIMMED_JOBS"       SEN_FASTQC_THREADS_PER_JOB="$SEN_FASTQC_THREADS_PER_JOB"       SEN_COVERAGE_JOBS="$SEN_COVERAGE_JOBS"       bash "$script"
}

run_stage "1/5 Download paired reads" "$REPO_ROOT/workflow/01_reads/01_download_fastq.sh"
run_stage "2/5 Raw FastQC + MultiQC" "$REPO_ROOT/workflow/02_qc/02_fastqc_raw.sh"
run_stage "3/5 fastp trimming" "$REPO_ROOT/workflow/03_typing/09_trim_fastp.sh"
run_stage "4/5 Trimmed FastQC + MultiQC" "$REPO_ROOT/workflow/02_qc/03_fastqc_trimmed.sh"
run_stage "5/5 Coverage estimate" "$REPO_ROOT/workflow/02_qc/04_estimate_coverage.sh"

expected_samples="$(grep -cve '^[[:space:]]*$' "$SRR_LIST")"
expected_fastq=$((expected_samples * 2))
raw_count="$(find "$SEN_ROOT/fastq" -type f -name '*.fastq.gz' 2>/dev/null | wc -l)"
trim_count="$(find "$SEN_ROOT/trimmed_fastq" -type f -name '*_trimmed_*.fastq.gz' 2>/dev/null | wc -l)"

echo
echo "=================================================="
echo " SMOKE READ/QC SUMMARY"
echo "=================================================="
echo "Expected paired FASTQs: $expected_fastq"
echo "Raw FASTQs found:       $raw_count"
echo "Trimmed FASTQs found:   $trim_count"
echo "Coverage table:         $SEN_ROOT/trimmed_fastq_qc/all_samples_coverage.csv"

if [[ "$raw_count" -ne "$expected_fastq" ]]; then
  echo "[FAIL] Raw FASTQ count mismatch." >&2
  exit 1
fi

if [[ "$trim_count" -ne "$expected_fastq" ]]; then
  echo "[FAIL] Trimmed FASTQ count mismatch." >&2
  exit 1
fi

echo
cat "$SEN_ROOT/trimmed_fastq_qc/all_samples_coverage.csv"
echo
echo "SMOKE READ/QC RESULT: PASS"
