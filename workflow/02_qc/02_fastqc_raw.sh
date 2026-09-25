#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd fastqc
require_cmd multiqc
require_cmd parallel

echo "=== RAW QC STARTED: $(date) ==="

RAW="${SEN_RAW_READS}"
QC="${SEN_RAW_QC_DIR:-$SEN_ROOT/raw_fastq_qc}"
LOGS="${SEN_LOG_DIR:-$SEN_ROOT/logs}"

mkdir -p "$QC" "$LOGS"

mapfile -t FASTQ_FILES < <(find "$RAW" -type f -name "*.fastq.gz")

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No .fastq.gz files found in $RAW or its subdirectories."
    exit 1
fi

echo "[INFO] Found ${#FASTQ_FILES[@]} FASTQ files."

JOBS="${SEN_FASTQC_RAW_JOBS:-30}"
THREADS_PER_JOB="${SEN_FASTQC_THREADS_PER_JOB:-2}"
MEM_PER_JOB="${SEN_FASTQC_MEM_PER_JOB:-2048m}"

export _JAVA_OPTIONS="-Xmx${MEM_PER_JOB}"

parallel -j "$JOBS" --eta "
    fastqc {} --threads ${THREADS_PER_JOB} --outdir ${QC}
" ::: "${FASTQ_FILES[@]}"

multiqc "$QC" -o "$QC"

echo "=== RAW QC COMPLETE: $(date) ==="
