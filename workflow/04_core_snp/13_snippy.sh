#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "\${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd snippy
require_cmd snippy-core
require_cmd snpEff
require_cmd snippy-vcf_to_tab
require_cmd parallel

PURE_DIR="\${SEN_KRAKEN_CLEAN_DIR:-$SEN_ROOT/Kraken_cleanup/clean_trimmed_fastq}"
OUT_DIR="\${SEN_SNIPPY_OUT:-$SEN_ROOT/Snippy_output}"
REF_FNA="\${SEN_REFERENCE_FNA:-$SEN_REFERENCE_DIR/reference.fna}"
REF_GFF="\${SEN_REFERENCE_GFF:-$SEN_REFERENCE_DIR/reference.gff}"
DB_DIR="\${SEN_SNPEFF_DIR:-$SEN_ROOT/snpEff_manual}"
PASS_LIST="\${SEN_PASS_LIST:-$SEN_ROOT/confirmed_enteritidis_list.txt}"

require_dir "$PURE_DIR"
require_file "$REF_FNA"
require_file "$REF_GFF"
require_file "$PASS_LIST"

JOBS="\${SEN_SNIPPY_JOBS:-60}"
THREADS_PER_JOB="\${SEN_SNIPPY_THREADS_PER_JOB:-2}"

mkdir -p "$OUT_DIR"

run_smart_snippy() {
    local R1=$1 REF=$2 OUT=$3 CPUS=$4
    local SN R2 SAMPLE_DIR
    SN=$(basename "$R1" _pure_1.fastq.gz)
    R2="\${R1/_pure_1.fastq.gz/_pure_2.fastq.gz}"
    SAMPLE_DIR="\${OUT}/\${SN}"

    if [[ -f "\${SAMPLE_DIR}/snps.tab" ]]; then
        return 0
    fi

    snippy --cpus "$CPUS" --outdir "$SAMPLE_DIR" --ref "$REF" \
           --pe1 "$R1" --pe2 "$R2" --cleanup --quiet
}
export -f run_smart_snippy

mapfile -t VALID_SAMPLES < "$PASS_LIST"
R1_TO_RUN=()
for SN in "\${VALID_SAMPLES[@]}"; do
    FILE="\${PURE_DIR}/\${SN}_pure_1.fastq.gz"
    [[ -f "$FILE" ]] && R1_TO_RUN+=("$FILE")
done

echo "[1/4] Checking/resuming Snippy alignments for \${#R1_TO_RUN[@]} genomes..."
parallel -j "$JOBS" --bar run_smart_snippy {} "$REF_FNA" "$OUT_DIR" "$THREADS_PER_JOB" ::: "\${R1_TO_RUN[@]}"

echo "[2/4] Building filtered core alignment..."
cd "$OUT_DIR"

PASSED_DIRS=()
for SN in "\${VALID_SAMPLES[@]}"; do
    [[ -d "\${OUT_DIR}/\${SN}" ]] && PASSED_DIRS+=("\${OUT_DIR}/\${SN}")
done

rm -f senbio_core.*
snippy-core --prefix senbio_core --ref "$REF_FNA" "\${PASSED_DIRS[@]}"

echo "[3/4] Running functional annotation on core VCF..."
snpEff -c "\${DB_DIR}/snpeff.config" -dataDir "\${DB_DIR}/data" -v ref \
       senbio_core.vcf > senbio_core.annotated.vcf

echo "[4/4] Generating final master report..."
snippy-vcf_to_tab --vcf senbio_core.annotated.vcf \
                  --ref "$REF_FNA" \
                  --gff "$REF_GFF" > SENBIO_FINAL_STRINGENT_REPORT.tab

echo "PIPELINE COMPLETE: $(date)"
echo "Final Samples in Alignment: $(grep -c '^>' senbio_core.full.aln)"
