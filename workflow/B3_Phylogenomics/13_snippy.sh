#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd snippy
require_cmd snippy-core
require_cmd parallel

PURE_DIR="${SEN_KRAKEN_CLEAN_DIR:-$SEN_ROOT/Kraken_cleanup/clean_trimmed_fastq}"
OUT_DIR="${SEN_SNIPPY_OUT:-$SEN_ROOT/Snippy_output}"
REF_FNA="${SEN_REFERENCE_FNA:-$SEN_REFERENCE_DIR/reference.fna}"
PASS_LIST="${SEN_PASS_LIST:-$SEN_ROOT/confirmed_enteritidis_list.txt}"

require_dir "$PURE_DIR"
require_file "$REF_FNA"
require_file "$PASS_LIST"

JOBS="${SEN_SNIPPY_JOBS:-60}"
THREADS_PER_JOB="${SEN_SNIPPY_THREADS_PER_JOB:-2}"

mkdir -p "$OUT_DIR"

run_smart_snippy() {
    local R1="$1" REF="$2" OUT="$3" CPUS="$4"
    local SN R2 SAMPLE_DIR
    SN="$(basename "$R1" _pure_1.fastq.gz)"
    R2="${R1/_pure_1.fastq.gz/_pure_2.fastq.gz}"
    SAMPLE_DIR="${OUT}/${SN}"

    if [[ -s "${SAMPLE_DIR}/snps.tab" && -s "${SAMPLE_DIR}/snps.vcf" ]]; then
        echo "[SKIP] Snippy $SN"
        return 0
    fi

    [[ -s "$R2" ]] || { echo "[ERROR] Missing R2 for $SN" >&2; return 1; }

    rm -rf "$SAMPLE_DIR"
    snippy       --cpus "$CPUS"       --outdir "$SAMPLE_DIR"       --ref "$REF"       --pe1 "$R1"       --pe2 "$R2"       --mincov 10       --minqual 100       --mapqual 60       --basequal 13       --minfrac 0       --cleanup       --quiet
    [[ -s "${SAMPLE_DIR}/snps.tab" && -s "${SAMPLE_DIR}/snps.vcf" ]]
}
export -f run_smart_snippy

mapfile -t VALID_SAMPLES < <(grep -vE '^[[:space:]]*(#|$)' "$PASS_LIST" | tr -d '\r')
R1_TO_RUN=()
for SN in "${VALID_SAMPLES[@]}"; do
    FILE="${PURE_DIR}/${SN}_pure_1.fastq.gz"
    [[ -s "$FILE" ]] || { echo "[ERROR] Missing cleaned R1 for $SN: $FILE" >&2; exit 1; }
    R1_TO_RUN+=("$FILE")
done

echo "[1/2] Checking/resuming Snippy alignments for ${#R1_TO_RUN[@]} genomes..."
parallel -j "$JOBS" --bar   run_smart_snippy {} "$REF_FNA" "$OUT_DIR" "$THREADS_PER_JOB"   ::: "${R1_TO_RUN[@]}"

echo "[2/2] Building core alignment..."
PASSED_DIRS=()
for SN in "${VALID_SAMPLES[@]}"; do
    [[ -s "${OUT_DIR}/${SN}/snps.tab" ]] || {
      echo "[ERROR] Missing completed Snippy output for $SN" >&2
      exit 1
    }
    PASSED_DIRS+=("${OUT_DIR}/${SN}")
done

cd "$OUT_DIR"
rm -f senbio_core.*
snippy-core --prefix senbio_core --ref "$REF_FNA" "${PASSED_DIRS[@]}"

[[ -s "$OUT_DIR/senbio_core.full.aln" && -s "$OUT_DIR/senbio_core.vcf" ]] || {
  echo "[ERROR] snippy-core did not create expected outputs." >&2
  exit 1
}

echo "SNIPPY CORE COMPLETE: $(date)"
echo "Alignment sequences: $(grep -c '^>' "$OUT_DIR/senbio_core.full.aln" || true)"
