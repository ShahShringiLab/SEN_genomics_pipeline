#!/bin/bash
set -euo pipefail

# 1. SET PATHS
BASE_DIR="/home/samuelajulo/SENBio/Final"
INPUT_DIR="$BASE_DIR/trimmed_fastq"
FAIL_DIR="$BASE_DIR/qc-failed_fastq"
COV_FILE="$BASE_DIR/trimmed_fastq_qc/all_samples_coverage.csv"
SISTR_FILE="$BASE_DIR/sistr_results_run/sistr_master_summary.csv"
META_FILE="$BASE_DIR/SEN_Genomes.csv"

mkdir -p "$FAIL_DIR"

echo "[STEP 1] Filtering by Coverage (>= 30X)..."
# Extract IDs from Coverage CSV (standard format)
COV_COL=$(head -n 1 "$COV_FILE" | tr ',' '\n' | grep -nx "FoldCoverage" | cut -d: -f1)
awk -F',' -v col="$COV_COL" 'NR > 1 && $col >= 30 {id=$1; gsub(/_trimmed(_1|_2)?/, "", id); gsub(/\.fastq\.gz/, "", id); print id}' "$COV_FILE" | sort -u > cov_pass_ids.tmp

# Get list of current files
ls "$INPUT_DIR"/*_trimmed_1.fastq.gz 2>/dev/null | xargs -n 1 basename | sed 's/_trimmed_1.fastq.gz//' | sort -u > current_before_cov.tmp
grep -vFf cov_pass_ids.tmp current_before_cov.tmp > cov_fail_list.tmp || true

# Move Coverage Failures
while read -r srr; do
    mv "$INPUT_DIR/${srr}_trimmed_1.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
    mv "$INPUT_DIR/${srr}_trimmed_2.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
done < cov_fail_list.tmp

echo "[STEP 2] Filtering by SISTR (Enteritidis)..."
# THE ROBUST FIX: Match "Enteritidis" and "SRR" pattern on the same line
grep "Enteritidis" "$SISTR_FILE" | grep -o "SRR[0-9]\+" | sort -u > sistr_pass_ids.tmp

# Check remaining files
ls "$INPUT_DIR"/*_trimmed_1.fastq.gz 2>/dev/null | xargs -n 1 basename | sed 's/_trimmed_1.fastq.gz//' | sort -u > current_after_cov.tmp
grep -vFf sistr_pass_ids.tmp current_after_cov.tmp > sistr_fail_list.tmp || true

# Move SISTR Failures
while read -r srr; do
    mv "$INPUT_DIR/${srr}_trimmed_1.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
    mv "$INPUT_DIR/${srr}_trimmed_2.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
done < sistr_fail_list.tmp

echo "[STEP 3] Updating Metadata..."
ls "$INPUT_DIR"/*_trimmed_1.fastq.gz 2>/dev/null | xargs -n 1 basename | sed 's/_trimmed_1.fastq.gz//' | sort -u > final_passed_list.tmp
# Sync CSV
HEADER=$(head -n 1 "$META_FILE")
RUN_COL=$(echo "$HEADER" | tr ',' '\n' | grep -nx "Run" | cut -d: -f1)
echo "$HEADER" > "$META_FILE.new"
awk -F',' -v col="$RUN_COL" 'NR==FNR{a[$1];next} $col in a' final_passed_list.tmp "$META_FILE" >> "$META_FILE.new"
mv "$META_FILE.new" "$META_FILE"

echo "DONE. Remaining: $(wc -l < final_passed_list.tmp)"
rm *.tmp
