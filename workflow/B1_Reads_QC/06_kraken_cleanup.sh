#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd kraken2
require_cmd seqkit
require_cmd parallel
require_cmd pigz
require_cmd python3

DB_PATH="${SEN_KRAKEN_DB:-}"
[[ -n "$DB_PATH" ]] || { echo "[ERROR] SEN_KRAKEN_DB is not set." >&2; exit 1; }
require_dir "$DB_PATH"

INPUT_DIR="${SEN_TRIMMED_READS}"
BASE_OUT="${SEN_KRAKEN_OUT:-$SEN_ROOT/Kraken_cleanup}"
CLEAN_DIR="$BASE_OUT/clean_trimmed_fastq"
REPORT_DIR="$BASE_OUT/kraken_reports"
FINAL_LOG="$BASE_OUT/full_db_results.csv"

PYTHON_EXEC="${SEN_PYTHON_EXEC:-python3}"
EXTRACT_SCRIPT="${SEN_KRAKEN_EXTRACT_SCRIPT:-$(command -v extract_kraken_reads.py || true)}"

if [[ -z "$EXTRACT_SCRIPT" ]]; then
    echo "[ERROR] extract_kraken_reads.py not found in PATH." >&2
    echo "[ERROR] Install/use environments/02_kraken.yaml, which includes KrakenTools." >&2
    exit 1
fi

GENOME_SIZE="${SEN_GENOME_SIZE_BP:-4800000}"
PARALLEL_JOBS="${SEN_KRAKEN_PARALLEL_JOBS:-24}"
KRAKEN_THREADS="${SEN_KRAKEN_THREADS:-4}"
PIGZ_THREADS="${SEN_PIGZ_THREADS:-2}"
TARGET_TAXID="${SEN_KRAKEN_TARGET_TAXID:-543}"
TMP_KRAKEN_DIR="${SEN_KRAKEN_TMP:-/dev/shm}"

mkdir -p "$CLEAN_DIR" "$REPORT_DIR" "$BASE_OUT"
echo "Sample,Initial_Depth,Final_Depth,Depth_Loss,Purity_Percent" > "$FINAL_LOG"

extract_target_percent() {
    awk -F'\t' -v tid="$TARGET_TAXID" '$5 == tid {print $1; exit}' "$1"
}
get_total_bases() {
    seqkit stats -T "$1" 2>/dev/null | awk 'NR==2 {print $5}'
}
export -f extract_target_percent get_total_bases

process_sample() {
    local fastq1="$1"
    local srr
    srr=$(basename "$fastq1" _trimmed_1.fastq.gz)
    local fastq2="${INPUT_DIR}/${srr}_trimmed_2.fastq.gz"
    local report_file="${REPORT_DIR}/${srr}.report"
    local clean1="${CLEAN_DIR}/${srr}_pure_1.fastq"
    local clean2="${CLEAN_DIR}/${srr}_pure_2.fastq"
    local kraken_out="${TMP_KRAKEN_DIR}/${srr}.kraken.$$.out"

    local b1 b2 initial_bases initial_depth
    b1=$(get_total_bases "$fastq1")
    b2=$(get_total_bases "$fastq2")
    initial_bases=$((b1 + b2))
    initial_depth=$(awk -v b="$initial_bases" -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f", b/g}')

    kraken2 --db "$DB_PATH" --threads "$KRAKEN_THREADS" --paired --gzip-compressed \
      --memory-mapping --output "$kraken_out" --report "$report_file" "$fastq1" "$fastq2"

    "$PYTHON_EXEC" "$EXTRACT_SCRIPT" -k "$kraken_out" -r "$report_file" \
      -s1 "$fastq1" -s2 "$fastq2" -t "$TARGET_TAXID" --include-children --fastq-output \
      -o "$clean1" -o2 "$clean2"

    rm -f "$kraken_out"

    if [[ -s "$clean1" ]]; then
        pigz -f -p "$PIGZ_THREADS" "$clean1" "$clean2"
        local f1 f2 final_bases final_depth loss purity
        f1=$(get_total_bases "${clean1}.gz")
        f2=$(get_total_bases "${clean2}.gz")
        final_bases=$((f1 + f2))
        final_depth=$(awk -v b="$final_bases" -v g="$GENOME_SIZE" 'BEGIN{printf "%.2f", b/g}')
        loss=$(awk -v i="$initial_depth" -v f="$final_depth" 'BEGIN{printf "%.2f", i-f}')
        purity=$(extract_target_percent "$report_file")
        purity=${purity:-0}
        echo "$srr,$initial_depth,$final_depth,$loss,$purity" >> "$FINAL_LOG"
    else
        echo "$srr,$initial_depth,0,$initial_depth,0" >> "$FINAL_LOG"
    fi
}

export -f process_sample
export DB_PATH INPUT_DIR REPORT_DIR CLEAN_DIR FINAL_LOG GENOME_SIZE KRAKEN_THREADS PIGZ_THREADS PYTHON_EXEC EXTRACT_SCRIPT TARGET_TAXID TMP_KRAKEN_DIR

find "$INPUT_DIR" -maxdepth 1 -name "*_trimmed_1.fastq.gz" -print0 | \
  parallel -0 -j "$PARALLEL_JOBS" --progress process_sample {}

generate_stats() {
    local col=$1 label=$2
    awk -F',' -v c="$col" 'NR>1 {
        if($c >= 30) a++; else if($c >= 25) b++; else if($c >= 20) d++; else e++
    } END {
        printf "%-15s | >30x: %d | 25-29x: %d | 20-24x: %d | <20x: %d\n", "'"$label"'", a, b, d, e
    }' "$FINAL_LOG"
}

generate_stats 2 "PRE-CLEANING"
generate_stats 3 "POST-CLEANING"
echo "[+] Done. Results saved to $FINAL_LOG"
