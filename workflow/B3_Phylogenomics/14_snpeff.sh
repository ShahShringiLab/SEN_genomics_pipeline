#!/usr/bin/env bash
# ==============================================================================
# SNPEFF DATABASE BUILDER - MINIMAL CONFIG VERSION
# ==============================================================================
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd snpEff
require_cmd snippy-vcf_to_tab
require_cmd parallel

# 1. CONFIGURATION
DB_DIR="${SEN_SNPEFF_DIR:-$SEN_ROOT/snpEff_manual}"
OUT_DIR="${SEN_SNIPPY_OUT:-$SEN_ROOT/Snippy_output}"
REF_FA="${SEN_REFERENCE_FA:-${SEN_REFERENCE_FNA:-$SEN_REFERENCE_DIR/reference.fna}}"
REF_GBK="${SEN_REFERENCE_GBK:-$SEN_REFERENCE_DIR/reference.gbk}"
REF_GFF="${SEN_REFERENCE_GFF:-$SEN_REFERENCE_DIR/reference.gff}"

require_file "$REF_FA"
require_file "$REF_GBK"

# 2. STEP 1: BUILD THE DATABASE
echo "[1/2] Building SnpEff database from your .gbk reference..."
mkdir -p "${DB_DIR}/data/ref"

# Copy files to SnpEff internal structure
cp "$REF_GBK" "${DB_DIR}/data/ref/genes.gbk"
cp "$REF_FA" "${DB_DIR}/data/ref/sequences.fa"

# Create a minimal config. We define ONLY the data directory and the genome name.
# We are NOT defining the codon table here to let SnpEff auto-detect from the GBK.
cat << EOF > "${DB_DIR}/snpeff.config"
data.dir = ${DB_DIR}/data
ref.genome : ref
EOF

# Build the database
# -v is verbose: it will tell us if it finds the genes.
snpEff build -c "${DB_DIR}/snpeff.config" -genbank -v ref

# 3. STEP 2: PARALLEL ANNOTATION FUNCTION
annotate_sample() {
    local SAMPLE_DIR=$1
    local DB_CONFIG=$2
    local REF_FA=$3
    local REF_GFF=$4

    local RAW_VCF="${SAMPLE_DIR}/snps.vcf"
    local ANN_VCF="${SAMPLE_DIR}/snps.annotated.vcf"
    local ANN_TAB="${SAMPLE_DIR}/snps.annotated.tab"

    if [ -f "$RAW_VCF" ]; then
        # Run SnpEff annotation
        snpEff -Xmx2g -c "$DB_CONFIG" -v ref "$RAW_VCF" > "$ANN_VCF"
        
        # Optional annotated TAB export. The publication SNP parser consumes
        # the annotated VCF directly, so a GFF is not required for a fresh run.
        if [[ -s "$REF_GFF" ]]; then
            snippy-vcf_to_tab --vcf "$ANN_VCF" --ref "$REF_FA" --gff "$REF_GFF" > "$ANN_TAB"
        else
            echo "[INFO] Reference GFF not present; skipping optional $ANN_TAB" >&2
        fi
    fi
}
export -f annotate_sample

echo "[2/2] Starting parallel annotation..."
find "$OUT_DIR" -maxdepth 1 -mindepth 1 -type d | \
parallel -j 16 --bar annotate_sample {} "${DB_DIR}/snpeff.config" "$REF_FA" "$REF_GFF"

echo "======================================================================"
echo "SUCCESS: All samples annotated."
echo "======================================================================"
