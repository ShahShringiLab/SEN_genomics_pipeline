#!/bin/bash
set -euo pipefail

echo "=== FASTQ DOWNLOAD STARTED: $(date) ==="
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate senbio

# directories
WD=/home/samuelajulo/SENBio/Final
SRR_LIST=$WD/SRR_list.clean.txt
SRA_DIR=$WD/sra_raw
FASTQ_DIR=$WD/fastq
TMP_DIR=$WD/tmp
LOGS=$WD/logs

mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMP_DIR" "$LOGS"

# Tune for your Threadripper (balanced CPU+IO)
PREFETCH_JOBS=24          # network bound
DUMP_JOBS=12              # CPU+IO bound
THREADS_PER_DUMP=12       # 12*12=144 threads total (good for 128C)

echo "[INFO] Prefetch jobs: $PREFETCH_JOBS"
echo "[INFO] fasterq-dump jobs: $DUMP_JOBS"
echo "[INFO] Threads per dump: $THREADS_PER_DUMP"

###########################################
# STEP 1 — PARALLEL PREFETCH (.sra files)
###########################################
echo "[STEP 1] Downloading SRA files in parallel..."

cat "$SRR_LIST" | parallel -j "$PREFETCH_JOBS" --eta --linebuffer "
  echo '[DL] {}'
  prefetch {} --output-directory '$SRA_DIR' >> '$LOGS/prefetch.log' 2>&1 || echo '[WARN] prefetch failed {}' >> '$LOGS/prefetch_failed.txt'
"

###########################################
# STEP 2 — PARALLEL CONVERSION TO FASTQ
###########################################
echo "[STEP 2] Extracting FASTQ using fasterq-dump..."

find "$SRA_DIR" -type f -name "*.sra" | parallel -j "$DUMP_JOBS" --eta --linebuffer "
  SRA={}
  SRR=\$(basename \"\$SRA\" .sra)
  OUTDIR='$FASTQ_DIR/'\"\$SRR\"

  # skip if already done
  if ls \"\$OUTDIR\"/*_1.fastq.gz >/dev/null 2>&1 || ls \"\$OUTDIR\"/*_1.fastq >/dev/null 2>&1; then
    echo '[SKIP] '\$SRR
    exit 0
  fi

  mkdir -p \"\$OUTDIR\"

  fasterq-dump \"\$SRA\" -O \"\$OUTDIR\" \
    --threads '$THREADS_PER_DUMP' \
    --temp '$TMP_DIR' \
    >> '$LOGS/fasterq.log' 2>&1 || { echo '[FAIL] fasterq-dump '\$SRR >> '$LOGS/fasterq_failed.txt'; exit 1; }

  pigz -p '$THREADS_PER_DUMP' -f \"\$OUTDIR\"/*.fastq
"

echo "=== FASTQ DOWNLOAD COMPLETE: $(date) ==="
echo "[INFO] Done. FASTQs in: $FASTQ_DIR"

