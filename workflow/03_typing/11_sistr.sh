#!/bin/bash
set -euo pipefail
# ==========================================================
# SISTR PIPELINE: cleaned_trimmed_fastq version (Kraken-cleaned)
# + CLEAN output folder before run
# ==========================================================

source ~/miniforge3/etc/profile.d/conda.sh
conda activate sistr_env

BASE="/home/samuelajulo/SENBio/Final"

# INPUT (Kraken-cleaned)
IN="$BASE/Kraken_cleanup/clean_trimmed_fastq"

OUT_ROOT="$BASE/sistr_results_run"
OUT_FASTA="$OUT_ROOT/mini_assemblies"
OUT_CSV="$OUT_ROOT/individual_csvs"
OUT_LOGS="$OUT_ROOT/skesa_logs"

echo "=== SISTR PIPELINE STARTED: $(date) ==="
echo "[INFO] IN: $IN"
echo "[INFO] OUT_ROOT: $OUT_ROOT"

# ----------------------------------------------------------
# CLEAN OUTPUT ROOT (fresh run)
# ----------------------------------------------------------
if [[ -d "$OUT_ROOT" ]]; then
  echo "[INFO] Removing existing output folder: $OUT_ROOT"
  rm -rf "$OUT_ROOT"
fi

mkdir -p "$OUT_FASTA" "$OUT_CSV" "$OUT_LOGS"

# Move into logs folder so per-sample SKESA logs land here
cd "$OUT_LOGS"

# Expect pairs like: SRRxxxxx_pure_1.fastq.gz and SRRxxxxx_pure_2.fastq.gz
find "$IN" -type f -name "*_pure_1.fastq.gz" -print0 | \
parallel -0 -j 20 --joblog "$OUT_ROOT/parallel_progress.log" '
  R1="{}"
  R2="${R1/_pure_1.fastq.gz/_pure_2.fastq.gz}"
  SN="$(basename "$R1" _pure_1.fastq.gz)"

  # Skip if pair missing
  if [[ ! -s "$R2" ]]; then
    echo "[SKIP] Missing R2 for $SN"
    exit 0
  fi

  # Step A: SKESA Assembly
  skesa --gz --fastq "$R1","$R2" --cores 4 --memory 8 \
    > "'"$OUT_FASTA"'/$SN.fasta" \
    2> "$SN.skesa.log"

  # Step B: SISTR Typing
  if [[ -s "'"$OUT_FASTA"'/$SN.fasta" ]]; then
    sistr -f csv -o "'"$OUT_CSV"'/$SN.csv" -n "$SN" "'"$OUT_FASTA"'/$SN.fasta" >/dev/null 2>&1
    echo "[DONE] $SN"
  else
    echo "[ERROR] Assembly failed for $SN"
    rm -f "'"$OUT_FASTA"'/$SN.fasta" || true
    exit 1
  fi
'

# Return to BASE to finish the report
cd "$BASE"

echo "=== GENERATING MASTER REPORT ==="
MASTER_REPORT="$OUT_ROOT/sistr_master_summary.csv"
FIRST_FILE="$(find "$OUT_CSV" -type f -name "*.csv" | head -n 1 || true)"

if [[ -n "${FIRST_FILE:-}" ]]; then
  head -n 1 "$FIRST_FILE" > "$MASTER_REPORT"
  find "$OUT_CSV" -type f -name "*.csv" -exec tail -n +2 {} + >> "$MASTER_REPORT"
  echo "SUCCESS: Master report created at $MASTER_REPORT"
else
  echo "ERROR: No SISTR results found."
fi

echo "=== PIPELINE FINISHED: $(date) ==="

