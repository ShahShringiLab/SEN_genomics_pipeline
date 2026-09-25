#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "\${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd skesa
require_cmd sistr
require_cmd parallel

BASE="$SEN_ROOT"
IN="\${SEN_KRAKEN_CLEAN_DIR:-$SEN_ROOT/Kraken_cleanup/clean_trimmed_fastq}"

OUT_ROOT="\${SEN_SISTR_OUT:-$SEN_ROOT/sistr_results_run}"
OUT_FASTA="$OUT_ROOT/mini_assemblies"
OUT_CSV="$OUT_ROOT/individual_csvs"
OUT_LOGS="$OUT_ROOT/skesa_logs"
JOBS="\${SEN_SISTR_JOBS:-20}"
SKESA_CORES="\${SEN_SKESA_CORES:-4}"
SKESA_MEMORY="\${SEN_SKESA_MEMORY_GB:-8}"

echo "=== SISTR PIPELINE STARTED: $(date) ==="
echo "[INFO] IN: $IN"
echo "[INFO] OUT_ROOT: $OUT_ROOT"

if [[ -d "$OUT_ROOT" ]]; then
  echo "[INFO] Removing existing output folder: $OUT_ROOT"
  rm -rf "$OUT_ROOT"
fi

mkdir -p "$OUT_FASTA" "$OUT_CSV" "$OUT_LOGS"
cd "$OUT_LOGS"

find "$IN" -type f -name "*_pure_1.fastq.gz" -print0 | parallel -0 -j "$JOBS" --joblog "$OUT_ROOT/parallel_progress.log" '
  R1="{}"
  R2="${R1/_pure_1.fastq.gz/_pure_2.fastq.gz}"
  SN="$(basename "$R1" _pure_1.fastq.gz)"

  if [[ ! -s "$R2" ]]; then
    echo "[SKIP] Missing R2 for $SN"
    exit 0
  fi

  skesa --gz --fastq "$R1","$R2" --cores '"$SKESA_CORES"' --memory '"$SKESA_MEMORY"'     > "'"$OUT_FASTA"'/$SN.fasta"     2> "$SN.skesa.log"

  if [[ -s "'"$OUT_FASTA"'/$SN.fasta" ]]; then
    sistr -f csv -o "'"$OUT_CSV"'/$SN.csv" -n "$SN" "'"$OUT_FASTA"'/$SN.fasta" >/dev/null 2>&1
    echo "[DONE] $SN"
  else
    echo "[ERROR] Assembly failed for $SN"
    rm -f "'"$OUT_FASTA"'/$SN.fasta" || true
    exit 1
  fi
'

cd "$BASE"

MASTER_REPORT="$OUT_ROOT/sistr_master_summary.csv"
FIRST_FILE="$(find "$OUT_CSV" -type f -name "*.csv" | head -n 1 || true)"

if [[ -n "\${FIRST_FILE:-}" ]]; then
  head -n 1 "$FIRST_FILE" > "$MASTER_REPORT"
  find "$OUT_CSV" -type f -name "*.csv" -exec tail -n +2 {} + >> "$MASTER_REPORT"
  echo "SUCCESS: Master report created at $MASTER_REPORT"
else
  echo "ERROR: No SISTR results found."
  exit 1
fi

echo "=== PIPELINE FINISHED: $(date) ==="
