#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd skesa
require_cmd sistr
require_cmd parallel

BASE="$SEN_ROOT"
IN="${SEN_KRAKEN_CLEAN_DIR:-$SEN_ROOT/Kraken_cleanup/clean_trimmed_fastq}"

OUT_ROOT="${SEN_SISTR_OUT:-$SEN_ROOT/sistr_results_run}"
OUT_FASTA="$OUT_ROOT/mini_assemblies"
OUT_CSV="$OUT_ROOT/individual_csvs"
OUT_LOGS="$OUT_ROOT/skesa_logs"
SISTR_LOGS="$OUT_ROOT/sistr_logs"
JOBS="${SEN_SISTR_JOBS:-20}"
SKESA_CORES="${SEN_SKESA_CORES:-4}"
SKESA_MEMORY="${SEN_SKESA_MEMORY_GB:-8}"
SISTR_THREADS="${SEN_SISTR_THREADS:-1}"

echo "=== SISTR PIPELINE STARTED: $(date) ==="
echo "[INFO] IN: $IN"
echo "[INFO] OUT_ROOT: $OUT_ROOT"

# Resume-safe: preserve valid per-sample outputs and recompute only incomplete samples.
mkdir -p "$OUT_FASTA" "$OUT_CSV" "$OUT_LOGS" "$SISTR_LOGS"

run_sample() {
  local R1="$1"
  local R2 SN FASTA CSV SKESA_LOG SISTR_LOG

  R2="${R1/_pure_1.fastq.gz/_pure_2.fastq.gz}"
  SN="$(basename "$R1" _pure_1.fastq.gz)"
  FASTA="$OUT_FASTA/$SN.fasta"
  CSV="$OUT_CSV/$SN.csv"
  SKESA_LOG="$OUT_LOGS/$SN.skesa.log"
  SISTR_LOG="$SISTR_LOGS/$SN.sistr.log"

  if [[ ! -s "$R2" ]]; then
    echo "[ERROR] Missing R2 for $SN" >&2
    return 1
  fi

  if [[ -s "$FASTA" ]]; then
    echo "[SKIP] SKESA $SN: assembly already present"
  else
    echo "[RUN] SKESA $SN"
    if ! skesa --gz --fastq "$R1","$R2"       --cores "$SKESA_CORES"       --memory "$SKESA_MEMORY"       > "$FASTA" 2> "$SKESA_LOG"; then
      rm -f "$FASTA"
      echo "[ERROR] SKESA failed for $SN; see $SKESA_LOG" >&2
      return 1
    fi
    [[ -s "$FASTA" ]] || {
      echo "[ERROR] SKESA produced an empty assembly for $SN" >&2
      return 1
    }
  fi

  if [[ -s "$CSV" ]]; then
    echo "[SKIP] SISTR $SN: result already present"
  else
    echo "[RUN] SISTR $SN"
    # Important: -n/--novel-alleles is NOT the genome-name option in SISTR.
    # Use -i <fasta> <genome_name> to assign the sample name explicitly.
    if ! sistr       -i "$FASTA" "$SN"       -f csv       -o "$CSV"       --qc       -t "$SISTR_THREADS"       > "$SISTR_LOG" 2>&1; then
      rm -f "$CSV"
      echo "[ERROR] SISTR failed for $SN; see $SISTR_LOG" >&2
      return 1
    fi
    [[ -s "$CSV" ]] || {
      echo "[ERROR] SISTR produced no CSV for $SN; see $SISTR_LOG" >&2
      return 1
    }
  fi

  echo "[DONE] $SN"
}

export -f run_sample
export OUT_FASTA OUT_CSV OUT_LOGS SISTR_LOGS SKESA_CORES SKESA_MEMORY SISTR_THREADS

find "$IN" -type f -name "*_pure_1.fastq.gz" -print0 |   parallel -0 -j "$JOBS" --joblog "$OUT_ROOT/parallel_progress.log" run_sample {}

MASTER_REPORT="$OUT_ROOT/sistr_master_summary.csv"
FIRST_FILE="$(find "$OUT_CSV" -type f -name "*.csv" -size +0c | head -n 1 || true)"

if [[ -n "${FIRST_FILE:-}" ]]; then
  head -n 1 "$FIRST_FILE" > "$MASTER_REPORT"
  find "$OUT_CSV" -type f -name "*.csv" -size +0c -exec tail -n +2 {} + >> "$MASTER_REPORT"
  echo "SUCCESS: Master report created at $MASTER_REPORT"
else
  echo "ERROR: No SISTR results found." >&2
  echo "[INFO] Inspect logs under: $SISTR_LOGS" >&2
  exit 1
fi

cd "$BASE"
echo "=== PIPELINE FINISHED: $(date) ==="
