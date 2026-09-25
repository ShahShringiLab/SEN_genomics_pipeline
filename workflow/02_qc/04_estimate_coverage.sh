#!/bin/bash

# 1. PATHS
INPUT_DIR="/home/samuelajulo/SENBio/Final/trimmed_fastq"
QC_DIR="/home/samuelajulo/SENBio/Final/trimmed_fastq_qc"
ALL_COV_FILE="$QC_DIR/all_samples_coverage.csv"
LOW_COV_FILE="$QC_DIR/low_coverage_samples_below_30x.csv"

mkdir -p "$QC_DIR"

# 2. CONSTANTS
GENOME_SIZE=4800000
THRESHOLD=30
CORES=64 

# 3. ENGINE FUNCTION (DYNAMIC READ LENGTH)
calc_cov() {
    local r1=$1
    local gs=$2
    
    local sample=$(basename "$r1" _trimmed_1.fastq.gz)
    local r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    
    # A. Calculate actual average read length from a sample of 10k reads
    # This handles 101, 150, 251, 301 bp automatically
    local avg_rl=$(zcat "$r1" | head -n 40000 | awk '{if(NR%4==2) {sum+=length($0); count++}} END {if(count>0) print int(sum/count); else print 0}')
    
    # B. Count total lines for exact total reads
    local line_count=$(zcat "$r1" "$r2" | wc -l)
    local total_reads=$((line_count / 4))
    
    # C. Accurate Math
    local total_bases=$(awk "BEGIN {print $total_reads * $avg_rl}")
    local coverage=$(awk "BEGIN {print $total_bases / $gs}")
    
    echo "$sample,$total_reads,$avg_rl,$total_bases,$coverage"
}
export -f calc_cov

# 4. EXECUTION
echo "--- Calculating Accurate Coverage (Dynamic Read Length) ---"

# Header for the CSV
echo "SampleID,ReadCount,AvgReadLen,TotalBases,FoldCoverage" > "$ALL_COV_FILE"

# Parallel processing
find "$INPUT_DIR" -name "*_1.fastq.gz" | \
parallel -j "$CORES" --line-buffer calc_cov {} "$GENOME_SIZE" >> "$ALL_COV_FILE"

# 5. FILTERING
echo "SampleID,FoldCoverage" > "$LOW_COV_FILE"
awk -F',' -v th="$THRESHOLD" 'NR > 1 && $5 < th {print $1 "," $5}' "$ALL_COV_FILE" >> "$LOW_COV_FILE"

echo -e "\n--- REPORT COMPLETE ---"
echo "Check $ALL_COV_FILE for dynamic read lengths per sample."
