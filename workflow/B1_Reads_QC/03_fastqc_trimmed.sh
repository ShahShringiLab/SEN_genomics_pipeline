#!/usr/bin/env bash
set -uo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd fastqc
require_cmd multiqc
require_cmd parallel

echo "=== TRIMMED QC (RESUME MODE) STARTED: $(date) ==="

TRIM_DIR="${SEN_TRIMMED_READS}"
QC_OUT="${SEN_TRIMMED_QC_DIR:-$SEN_ROOT/trimmed_fastq_qc}"
LOGS="${SEN_LOG_DIR:-$SEN_ROOT/logs}"

mkdir -p "$QC_OUT" "$LOGS"

mapfile -t ALL_FILES < <(find "$TRIM_DIR" -type f -name "*_trimmed_*.fastq.gz")
FILES_TO_RUN=()

echo "[INFO] Checking which of the ${#ALL_FILES[@]} files need processing..."

for f in "${ALL_FILES[@]}"; do
    BN=$(basename "$f")
    EXPECTED_OUT="${QC_OUT}/${BN/.fastq.gz/_fastqc.zip}"
    if [[ ! -f "$EXPECTED_OUT" ]]; then
        FILES_TO_RUN+=("$f")
    fi
done

if [[ ${#FILES_TO_RUN[@]} -eq 0 ]]; then
    echo "[INFO] All FastQC reports already exist. Skipping to MultiQC."
else
    JOBS="${SEN_FASTQC_TRIMMED_JOBS:-40}"
    THREADS_PER_JOB="${SEN_FASTQC_THREADS_PER_JOB:-2}"
    export _JAVA_OPTIONS="-Xmx${SEN_FASTQC_MEM_PER_JOB:-2048m}"

    parallel -j "$JOBS" --bar "
        fastqc {} --threads ${THREADS_PER_JOB} --outdir ${QC_OUT} || true
    " ::: "${FILES_TO_RUN[@]}"
fi

multiqc "$QC_OUT"     --outdir "$QC_OUT"     --filename "trimmed_data_multiqc_report"     --force     --verbose

echo "=== PROCESS COMPLETE: $(date) ==="
