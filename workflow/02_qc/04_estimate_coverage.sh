#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd parallel
require_cmd zcat
require_cmd awk

INPUT_DIR="${SEN_TRIMMED_READS}"
QC_DIR="${SEN_TRIMMED_QC_DIR:-$SEN_ROOT/trimmed_fastq_qc}"
ALL_COV_FILE="$QC_DIR/all_samples_coverage.csv"
LOW_COV_FILE="$QC_DIR/low_coverage_samples_below_30x.csv"

mkdir -p "$QC_DIR"

GENOME_SIZE="${SEN_GENOME_SIZE_BP:-4800000}"
THRESHOLD="${SEN_MIN_COVERAGE:-30}"
CORES="${SEN_COVERAGE_JOBS:-64}"

calc_cov() {
    local r1=$1
    local gs=$2

    local sample
    sample=$(basename "$r1" _trimmed_1.fastq.gz)
    local r2="${r1/_1.fastq.gz/_2.fastq.gz}"

    local avg_rl
    avg_rl=$(zcat "$r1" | head -n 40000 | awk '{if(NR%4==2) {sum+=length($0); count++}} END {if(count>0) print int(sum/count); else print 0}')

    local line_count
    line_count=$(zcat "$r1" "$r2" | wc -l)
    local total_reads=$((line_count / 4))

    local total_bases
    total_bases=$(awk "BEGIN {print $total_reads * $avg_rl}")
    local coverage
    coverage=$(awk "BEGIN {print $total_bases / $gs}")

    echo "$sample,$total_reads,$avg_rl,$total_bases,$coverage"
}
export -f calc_cov

echo "SampleID,ReadCount,AvgReadLen,TotalBases,FoldCoverage" > "$ALL_COV_FILE"

find "$INPUT_DIR" -name "*_1.fastq.gz" | parallel -j "$CORES" --line-buffer calc_cov {} "$GENOME_SIZE" >> "$ALL_COV_FILE"

echo "SampleID,FoldCoverage" > "$LOW_COV_FILE"
awk -F',' -v th="$THRESHOLD" 'NR > 1 && $5 < th {print $1 "," $5}' "$ALL_COV_FILE" >> "$LOW_COV_FILE"

echo "Check $ALL_COV_FILE for dynamic read lengths per sample."
