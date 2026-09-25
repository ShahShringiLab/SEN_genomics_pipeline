#!/bin/bash
set -euo pipefail

# ============================================================
# Kraken2 Cleanup with BEFORE vs AFTER Depth Distribution
# ============================================================

# -------------------------
# PATHS
# -------------------------
DB_PATH="/home/samuelajulo/kraken_db"
INPUT_DIR="/home/samuelajulo/SENBio/Final/trimmed_fastq"
BASE_OUT="/home/samuelajulo/SENBio/Final/Kraken_cleanup"

CLEAN_DIR="$BASE_OUT/clean_trimmed_fastq"
REPORT_DIR="$BASE_OUT/kraken_reports"
FINAL_LOG="$BASE_OUT/full_db_results.csv"

# -------------------------
# TOOL LOCATIONS
# -------------------------
PYTHON_EXEC="/home/samuelajulo/miniforge3/envs/contamination_qc/bin/python3"
EXTRACT_SCRIPT="/home/samuelajulo/SENBio/Final/extract_kraken_reads_custom.py"

# -------------------------
# SETTINGS
# -------------------------
GENOME_SIZE=4800000
PARALLEL_JOBS=24
KRAKEN_THREADS=4
PIGZ_THREADS=2
TARGET_TAXID="543"
TMP_KRAKEN_DIR="/dev/shm"

# -------------------------
# FRESH START
# -------------------------
echo "[!] Starting Fresh..."
mkdir -p "$CLEAN_DIR" "$REPORT_DIR" "$BASE_OUT"

# Header for the CSV
echo "Sample,Initial_Depth,Final_Depth,Depth_Loss,Purity_Percent" > "$FINAL_LOG"

# -------------------------
# EXPORTABLE HELPERS
# -------------------------
extract_target_percent() {
    # Extracts the % of reads assigned to Salmonella (TaxID 543)
    awk -F'\t' -v tid="$TARGET_TAXID" '$5 == tid {print $1; exit}' "$1"
}

get_total_bases() {
    # Uses seqkit to get total base count safely
    seqkit stats -T "$1" 2>/dev/null | awk 'NR==2 {print $5}'
}

export -f extract_target_percent
export -f get_total_bases

# -------------------------
# PROCESSOR FUNCTION
# -------------------------
process_sample() {
    local fastq1="$1"
    local srr=$(basename "$fastq1" _trimmed_1.fastq.gz)
    local fastq2="${INPUT_DIR}/${srr}_trimmed_2.fastq.gz"
    local report_file="${REPORT_DIR}/${srr}.report"
    local clean1="${CLEAN_DIR}/${srr}_pure_1.fastq"
    local clean2="${CLEAN_DIR}/${srr}_pure_2.fastq"
    local kraken_out="${TMP_KRAKEN_DIR}/${srr}.kraken.$$.out"

    # 1. Get Initial Depth
    local b1=$(get_total_bases "$fastq1")
    local b2=$(get_total_bases "$fastq2")
    local initial_bases=$((b1 + b2))
    local initial_depth=$(awk -v b="$initial_bases" -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f", b/g}')

    # 2. Run Kraken2
    kraken2 --db "$DB_PATH" --threads "$KRAKEN_THREADS" --paired --gzip-compressed \
            --memory-mapping --output "$kraken_out" --report "$report_file" \
            "$fastq1" "$fastq2"

    # 3. Extract Target Reads
    "$PYTHON_EXEC" "$EXTRACT_SCRIPT" -k "$kraken_out" --report "$report_file" \
        -s1 "$fastq1" -s2 "$fastq2" -t "$TARGET_TAXID" --include-children --fastq-output \
        -o "$clean1" -o2 "$clean2"
    
    rm -f "$kraken_out"

    # 4. Process Clean Files
    if [[ -s "$clean1" ]]; then
        pigz -f -p "$PIGZ_THREADS" "$clean1" "$clean2"
        
        local f1=$(get_total_bases "${clean1}.gz")
        local f2=$(get_total_bases "${clean2}.gz")
        local final_bases=$((f1 + f2))
        local final_depth=$(awk -v b="$final_bases" -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f", b/g}')
        
        local loss=$(awk -v i="$initial_depth" -v f="$final_depth" 'BEGIN{printf "%.2f", i-f}')
        local purity=$(extract_target_percent "$report_file")
        purity=${purity:-0}

        echo "$srr,$initial_depth,$final_depth,$loss,$purity" >> "$FINAL_LOG"
    else
        echo "$srr,$initial_depth,0,$initial_depth,0" >> "$FINAL_LOG"
    fi
}

export -f process_sample
export DB_PATH INPUT_DIR REPORT_DIR CLEAN_DIR FINAL_LOG GENOME_SIZE KRAKEN_THREADS PIGZ_THREADS PYTHON_EXEC EXTRACT_SCRIPT TARGET_TAXID TMP_KRAKEN_DIR

# -------------------------
# EXECUTION
# -------------------------
echo "[>] Running Kraken2 Cleanup on Parallel Jobs..."
find "$INPUT_DIR" -maxdepth 1 -name "*_trimmed_1.fastq.gz" -print0 | \
    parallel -0 -j "$PARALLEL_JOBS" --progress process_sample {}

# -------------------------
# DISTRIBUTION SUMMARY
# -------------------------
echo -e "\n-------------------------------------------------------"
echo "  DEPTH DISTRIBUTION SUMMARY (BEFORE VS AFTER)"
echo "-------------------------------------------------------"

generate_stats() {
    local col=$1
    local label=$2
    awk -F',' -v c="$col" 'NR>1 {
        if($c >= 30) a++; 
        else if($c >= 25) b++; 
        else if($c >= 20) d++;
        else e++
    } END {
        printf "%-15s | >30x: %d | 25-29x: %d | 20-24x: %d | <20x: %d\n", "'"$label"'", a, b, d, e
    }' "$FINAL_LOG"
}

generate_stats 2 "PRE-CLEANING"
generate_stats 3 "POST-CLEANING"

echo "-------------------------------------------------------"
echo "[+] Done. Results saved to $FINAL_LOG"
