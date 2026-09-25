#!/bin/bash
set -uo pipefail

echo "=== TRIMMED QC (RESUME MODE) STARTED: $(date) ==="

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate senbio

# =============================
# PATHS
# =============================
WD=/home/samuelajulo/SENBio/Final
TRIM_DIR=${WD}/trimmed_fastq
QC_OUT=${WD}/trimmed_fastq_qc
LOGS=${WD}/logs

mkdir -p "$QC_OUT" "$LOGS"

# =============================
# FILE DISCOVERY
# =============================
ALL_FILES=($(find "$TRIM_DIR" -type f -name "*_trimmed_*.fastq.gz"))
FILES_TO_RUN=()

echo "[INFO] Checking which of the ${#ALL_FILES[@]} files need processing..."

for f in "${ALL_FILES[@]}"; do
    # Get just the filename (e.g., SRR123_trimmed_1.fastq.gz)
    BN=$(basename "$f")
    # FastQC output name is usually the filename minus .fastq.gz plus _fastqc.html
    # We check for the .zip file as it's the last thing FastQC creates
    EXPECTED_OUT="${QC_OUT}/${BN/.fastq.gz/_fastqc.zip}"

    if [ ! -f "$EXPECTED_OUT" ]; then
        FILES_TO_RUN+=("$f")
    fi
done

# =============================
# EXECUTION LOGIC
# =============================
if [ ${#FILES_TO_RUN[@]} -eq 0 ]; then
    echo "[INFO] All FastQC reports already exist. Skipping to MultiQC."
else
    echo "[INFO] ${#FILES_TO_RUN[@]} files still need QC. Running now..."
    
    JOBS=40
    THREADS_PER_JOB=2
    export _JAVA_OPTIONS="-Xmx2048m"

    parallel -j $JOBS --bar "
        fastqc {} --threads ${THREADS_PER_JOB} --outdir ${QC_OUT} || true
    " ::: "${FILES_TO_RUN[@]}"
fi

# =============================
# MULTIQC (ALWAYS RUNS)
# =============================
echo "[INFO] Running MultiQC..."

multiqc "$QC_OUT" \
    --outdir "$QC_OUT" \
    --filename "trimmed_data_multiqc_report" \
    --force \
    --verbose

echo "=== PROCESS COMPLETE: $(date) ==="
