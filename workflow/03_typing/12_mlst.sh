#!/bin/bash
# ==========================================================
# MLST PIPELINE: AUTO-INSTALL & HIGH-SPEED EXECUTION
# ==========================================================

source ~/miniforge3/etc/profile.d/conda.sh

# 1. Environment Setup
ENV_NAME="mlst_env"
if { conda env list | grep "$ENV_NAME"; } >/dev/null 2>&1; then
    echo "=== Activating existing $ENV_NAME ==="
    conda activate $ENV_NAME
else
    echo "=== Creating $ENV_NAME and installing mlst ==="
    conda create -n $ENV_NAME -c bioconda -c conda-forge mlst parallel -y
    conda activate $ENV_NAME
fi

# 2. Paths
BASE="/home/samuelajulo/SENBio/Final"
IN_FASTA="$BASE/sistr_results_run/mini_assemblies"
OUT_DIR="$BASE/mlst_results_run"
MASTER_REPORT="$OUT_DIR/mlst_master_report.csv"

mkdir -p "$OUT_DIR"

echo "=== MLST PIPELINE STARTED: $(date) ==="

# 3. System Optimization
# We detect available CPU cores to maximize parallelization
CORES=$(nproc)
echo "Utilizing $CORES cores for analysis..."

# 4. Execution Logic
# Batching files (-n 100) is key to speed for 3400+ files
echo "Scanning $IN_FASTA for assemblies..."

find "$IN_FASTA" -name "*.fasta" | parallel -j $CORES -n 100 "mlst --threads 1 --csv {}" > "$OUT_DIR/all_results.raw.csv"

# 5. Header Addition for R Compatibility
# Traditional Achtman 7-gene scheme headers
echo "file,species,ST,aroC,dnaN,hemD,hisD,purE,sucA,thrA" > "$MASTER_REPORT"
cat "$OUT_DIR/all_results.raw.csv" >> "$MASTER_REPORT"

# Cleanup
rm "$OUT_DIR/all_results.raw.csv"

echo "=== SUCCESS: Master MLST report created at $MASTER_REPORT ==="
echo "=== PIPELINE FINISHED: $(date) ==="
