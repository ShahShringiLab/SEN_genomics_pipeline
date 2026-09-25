#!/bin/bash
# ==============================================================================
# ULTIMATE SNIPPY PIPELINE: CLEANED & OPTIMIZED
# ==============================================================================
set -o pipefail
# set -u is disabled briefly here to prevent Conda activation errors
set +u 

# 1. ENVIRONMENT & PATHS
source ~/miniforge3/etc/profile.d/conda.sh
conda activate snippy_env
set -u # Re-enable strict variable checking

WD="/home/samuelajulo/SENBio/Final"
PURE_DIR="${WD}/Kraken_cleanup/clean_trimmed_fastq"
OUT_DIR="${WD}/Snippy_output"
REF_FNA="${WD}/reference/reference.fna"
REF_GBK="${WD}/reference/reference.gbk"
REF_GFF="${WD}/reference/reference.gff"
DB_DIR="${WD}/snpEff_manual"
PASS_LIST="${WD}/confirmed_enteritidis_list.txt"

# 2. SETTINGS
JOBS=60
THREADS_PER_JOB=2

# 3. DEFINE THE SMART ALIGNMENT FUNCTION
run_smart_snippy() {
    local R1=$1; local REF=$2; local OUT=$3; local CPUS=$4
    local SN=$(basename "$R1" _pure_1.fastq.gz)
    local R2="${R1/_pure_1.fastq.gz/_pure_2.fastq.gz}"
    local SAMPLE_DIR="${OUT}/${SN}"

    # SKIP LOGIC: If snps.tab exists, it's finished.
    if [ -f "${SAMPLE_DIR}/snps.tab" ]; then
        return 0
    else
        snippy --cpus "$CPUS" --outdir "$SAMPLE_DIR" --ref "$REF" \
               --pe1 "$R1" --pe2 "$R2" --cleanup --quiet
    fi
}
export -f run_smart_snippy

# 4. STEP 1: ALIGNMENT (ONLY FOR PASSED GENOMES)
echo "[1/4] Checking/Resuming alignments for 3,310 confirmed genomes..."
mapfile -t VALID_SAMPLES < "$PASS_LIST"
R1_TO_RUN=()
for SN in "${VALID_SAMPLES[@]}"; do
    FILE="${PURE_DIR}/${SN}_pure_1.fastq.gz"
    [ -f "$FILE" ] && R1_TO_RUN+=("$FILE")
done

parallel -j "$JOBS" --bar run_smart_snippy {} "$REF_FNA" "$OUT_DIR" "$THREADS_PER_JOB" ::: "${R1_TO_RUN[@]}"

# 5. STEP 2: FILTERED CORE MERGE
echo "[2/4] Merging high-quality alignment (Count: 3311 expected)..."
cd "$OUT_DIR"
PASSED_DIRS=""
for SN in "${VALID_SAMPLES[@]}"; do
    [ -d "${OUT_DIR}/${SN}" ] && PASSED_DIRS+="${OUT_DIR}/${SN} "
done

# Clear old files to ensure the new 3311 count is saved correctly
rm -f senbio_core.* snippy-core --prefix senbio_core --ref "$REF_FNA" $PASSED_DIRS

# 6. STEP 3: ANNOTATION (snpEff)
echo "[3/4] Running functional annotation on core VCF..."
snpEff -c "${DB_DIR}/snpeff.config" -dataDir "${DB_DIR}/data" -v ref \
       senbio_core.vcf > senbio_core.annotated.vcf

# 7. STEP 4: FINAL REPORT
echo "[4/4] Generating final master report..."
snippy-vcf_to_tab --vcf senbio_core.annotated.vcf \
                  --ref "$REF_FNA" \
                  --gff "$REF_GFF" > SENBIO_FINAL_STRINGENT_REPORT.tab

echo "======================================================================"
echo "PIPELINE COMPLETE: $(date)"
echo "Final Samples in Alignment: $(grep -c '>' senbio_core.full.aln)"
echo "======================================================================"
