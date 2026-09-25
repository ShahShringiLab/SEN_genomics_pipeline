#!/bin/bash
set -euo pipefail

echo "=== RAW QC STARTED: $(date) ==="

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate senbio

# directory setup
WD=/home/samuelajulo/SENBio/Final
RAW=${WD}/fastq
QC=${WD}/raw_fastq_qc
LOGS=${WD}/logs

mkdir -p "$QC" "$LOGS"

# =============================
# FILE DISCOVERY (FIXED)
# =============================
# Using 'find' instead of 'ls' to search recursively inside SRR subfolders
# -type f searches only for files
# -name filters for fastq.gz
FASTQ_FILES=($(find "$RAW" -type f -name "*.fastq.gz"))

# Safety check: Exit if no files are found
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No .fastq.gz files found in $RAW or its subdirectories."
    echo "       Please verify your directory structure."
    exit 1
fi

echo "[INFO] Found ${#FASTQ_FILES[@]} FASTQ files."

# =============================
# PARALLEL FASTQC CONFIGURATION
# =============================
# Note: Ensure your system has (30 * 2) = 60 Threads and (30 * 2GB) = 60GB RAM available.
JOBS=30           
THREADS_PER_JOB=2
MEM_PER_JOB=2048m 

export _JAVA_OPTIONS="-Xmx${MEM_PER_JOB}"

echo "[INFO] Running FastQC in parallel:"
echo "      Jobs: $JOBS"
echo "      Threads/job: $THREADS_PER_JOB"
echo "      Memory/job: $MEM_PER_JOB"

# Run GNU Parallel
# The {} is replaced by each file found in the FASTQ_FILES array
parallel -j $JOBS --eta "
    fastqc {} \
      --threads ${THREADS_PER_JOB} \
      --outdir ${QC}
" ::: "${FASTQ_FILES[@]}"

echo "[INFO] Running MultiQC..."
multiqc "$QC" -o "$QC"

echo "=== RAW QC COMPLETE: $(date) ==="
