#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

INPUT_DIR="${SEN_TRIMMED_READS}"
FAIL_DIR="${SEN_QC_FAIL_DIR:-$SEN_ROOT/qc-failed_fastq}"
COV_FILE="${SEN_COVERAGE_FILE:-$SEN_ROOT/trimmed_fastq_qc/all_samples_coverage.csv}"
SISTR_FILE="${SEN_SISTR_SUMMARY:-$SEN_ROOT/sistr_results_run/sistr_master_summary.csv}"
META_FILE="${SEN_METADATA_FILE:-$SEN_ROOT/metadata/SEN_Genomes.csv}"

require_file "$COV_FILE"
require_file "$SISTR_FILE"
require_file "$META_FILE"
mkdir -p "$FAIL_DIR"

echo "[STEP 1] Filtering by Coverage (>= 30X)..."
COV_COL=$(head -n 1 "$COV_FILE" | tr ',' '\n' | grep -nx "FoldCoverage" | cut -d: -f1)
awk -F',' -v col="$COV_COL" 'NR > 1 && $col >= 30 {id=$1; gsub(/_trimmed(_1|_2)?/, "", id); gsub(/\.fastq\.gz/, "", id); print id}' "$COV_FILE" | sort -u > cov_pass_ids.tmp

ls "$INPUT_DIR"/*_trimmed_1.fastq.gz 2>/dev/null | xargs -n 1 basename | sed 's/_trimmed_1.fastq.gz//' | sort -u > current_before_cov.tmp
grep -vFf cov_pass_ids.tmp current_before_cov.tmp > cov_fail_list.tmp || true

while read -r srr; do
    mv "$INPUT_DIR/${srr}_trimmed_1.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
    mv "$INPUT_DIR/${srr}_trimmed_2.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
done < cov_fail_list.tmp

echo "[STEP 2] Filtering by SISTR (Enteritidis)..."
grep "Enteritidis" "$SISTR_FILE" | grep -o "SRR[0-9]\+" | sort -u > sistr_pass_ids.tmp

ls "$INPUT_DIR"/*_trimmed_1.fastq.gz 2>/dev/null | xargs -n 1 basename | sed 's/_trimmed_1.fastq.gz//' | sort -u > current_after_cov.tmp
grep -vFf sistr_pass_ids.tmp current_after_cov.tmp > sistr_fail_list.tmp || true

while read -r srr; do
    mv "$INPUT_DIR/${srr}_trimmed_1.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
    mv "$INPUT_DIR/${srr}_trimmed_2.fastq.gz" "$FAIL_DIR/" 2>/dev/null || true
done < sistr_fail_list.tmp

echo "[STEP 3] Updating Metadata..."
ls "$INPUT_DIR"/*_trimmed_1.fastq.gz 2>/dev/null | xargs -n 1 basename | sed 's/_trimmed_1.fastq.gz//' | sort -u > final_passed_list.tmp

HEADER=$(head -n 1 "$META_FILE")
RUN_COL=$(echo "$HEADER" | tr ',' '\n' | grep -nx "Run" | cut -d: -f1)
echo "$HEADER" > "$META_FILE.new"
awk -F',' -v col="$RUN_COL" 'NR==FNR{a[$1];next} $col in a' final_passed_list.tmp "$META_FILE" >> "$META_FILE.new"
mv "$META_FILE.new" "$META_FILE"

echo "DONE. Remaining: $(wc -l < final_passed_list.tmp)"
rm -f ./*.tmp
