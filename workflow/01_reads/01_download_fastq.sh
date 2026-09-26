#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd prefetch
require_cmd fasterq-dump
require_cmd parallel
require_cmd pigz

echo "=== FASTQ DOWNLOAD STARTED: $(date) ==="

SRR_LIST="${SEN_SRR_LIST:-$SEN_ROOT/SRR_list.clean.txt}"
SRA_DIR="${SEN_SRA_DIR:-$SEN_ROOT/sra_raw}"
FASTQ_DIR="${SEN_RAW_READS}"
TMP_DIR="${SEN_TMP_DIR:-$SEN_ROOT/tmp}"
LOGS="${SEN_LOG_DIR:-$SEN_ROOT/logs}"

require_file "$SRR_LIST"
mkdir -p "$SRA_DIR" "$FASTQ_DIR" "$TMP_DIR" "$LOGS"

CPU_COUNT="${SEN_CPU_COUNT:-$(nproc 2>/dev/null || echo 4)}"
DEFAULT_PREFETCH_JOBS=$(( CPU_COUNT >= 16 ? 8 : (CPU_COUNT >= 8 ? 4 : 2) ))
DEFAULT_DUMP_JOBS=$(( CPU_COUNT >= 16 ? 4 : (CPU_COUNT >= 8 ? 2 : 1) ))
DEFAULT_THREADS_PER_DUMP=$(( CPU_COUNT / DEFAULT_DUMP_JOBS ))
(( DEFAULT_THREADS_PER_DUMP > 8 )) && DEFAULT_THREADS_PER_DUMP=8
(( DEFAULT_THREADS_PER_DUMP < 2 )) && DEFAULT_THREADS_PER_DUMP=2

PREFETCH_JOBS="${SEN_PREFETCH_JOBS:-$DEFAULT_PREFETCH_JOBS}"
DUMP_JOBS="${SEN_DUMP_JOBS:-$DEFAULT_DUMP_JOBS}"
THREADS_PER_DUMP="${SEN_THREADS_PER_DUMP:-$DEFAULT_THREADS_PER_DUMP}"

echo "[INFO] SEN_ROOT: $SEN_ROOT"
echo "[INFO] SRR list: $SRR_LIST"
echo "[INFO] Prefetch jobs: $PREFETCH_JOBS"
echo "[INFO] fasterq-dump jobs: $DUMP_JOBS"
echo "[INFO] Threads per dump: $THREADS_PER_DUMP"

cat "$SRR_LIST" | parallel -j "$PREFETCH_JOBS" --eta --linebuffer "
  echo '[DL] {}'
  prefetch {} --output-directory '$SRA_DIR' >> '$LOGS/prefetch.log' 2>&1 || echo '[WARN] prefetch failed {}' >> '$LOGS/prefetch_failed.txt'
"

find "$SRA_DIR" -type f -name "*.sra" | parallel -j "$DUMP_JOBS" --eta --linebuffer "
  SRA={}
  SRR=\$(basename \"\$SRA\" .sra)
  OUTDIR='$FASTQ_DIR/'\"\$SRR\"

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
