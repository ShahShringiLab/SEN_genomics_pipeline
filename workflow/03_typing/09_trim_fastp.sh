#!/bin/bash
set -euo pipefail

echo "=== OPTIMIZED TRIMMING STARTED: $(date) ==="

# =====================================
# 1. ENABLE CONDA
# =====================================
source ~/miniforge3/etc/profile.d/conda.sh
conda activate senbio

# =====================================
# 2. PATHS
# =====================================
WD=/home/samuelajulo/SENBio/Final
RAW=${WD}/fastq
TRIM=${WD}/trimmed_fastq
LOGS=${WD}/logs

mkdir -p "$TRIM" "$LOGS"

# =====================================
# 3. HARDWARE OPTIMIZATION (128 Cores)
# =====================================
# fastp scales best at 2 threads. 
# 60 jobs * 2 threads = 120 threads (leaves 8 for system/disk overhead)
JOBS=60
THREADS_PER_JOB=2

echo "[INFO] Parallel Jobs: $JOBS"
echo "[INFO] Threads per Job: $THREADS_PER_JOB"

# =====================================
# 4. AUTODETECT SAMPLES
# =====================================
R1_FILES=($(find "$RAW" -type f -name "*_1.fastq.gz"))

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No files ending in _1.fastq.gz found in $RAW"
    exit 1
fi

echo "[INFO] Found ${#R1_FILES[@]} samples to trim."

# =====================================
# 5. OPTIMIZED TRIMMING FUNCTION
# =====================================
trim_one_sample () {
    FULL_R1=$1
    TRIM_DIR=$2
    THREADS=$3
    LOG_DIR=$4

    FULL_R2="${FULL_R1/_1.fastq.gz/_2.fastq.gz}"
    BASENAME=$(basename "$FULL_R1" _1.fastq.gz)

    OUT1=${TRIM_DIR}/${BASENAME}_trimmed_1.fastq.gz
    OUT2=${TRIM_DIR}/${BASENAME}_trimmed_2.fastq.gz
    JSON=${LOG_DIR}/${BASENAME}.json
    HTML=${LOG_DIR}/${BASENAME}.html

    if [ ! -f "$FULL_R2" ]; then
        echo "[WARN] Read 2 not found for $BASENAME. Skipping."
        return
    fi

    # Improved fastp parameters:
    # --trim_poly_g: Removes G-tails common in NovaSeq/NextSeq
    # --cut_right: Sliding window trimming (scans 5' to 3')
    # --correction: PE base correction if reads overlap
    # -q 20: Higher quality threshold for kept bases
    fastp \
        -i "$FULL_R1" \
        -I "$FULL_R2" \
        -o "$OUT1" \
        -O "$OUT2" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --trim_poly_g --poly_g_min_len 10 \
        --cut_right --cut_window_size 4 --cut_mean_quality 20 \
        --correction \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --length_required 50 \
        --compression 4 \
        --json "$JSON" \
        --html "$HTML" \
        >/dev/null 2>&1
}
export -f trim_one_sample

# =====================================
# 6. RUN FASTP IN PARALLEL
# =====================================
echo "[INFO] Running optimized fastp trimming..."

parallel -j $JOBS --bar \
    trim_one_sample {} "$TRIM" "$THREADS_PER_JOB" "$LOGS" \
    ::: "${R1_FILES[@]}"

echo "=== TRIMMING COMPLETE: $(date) ==="
