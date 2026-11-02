#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd iqtree2

GUBBINS_DIR="${SEN_GUBBINS_OUT:-$SEN_ROOT/gubbin}"
INPUT_FILE="${SEN_IQTREE_INPUT:-$GUBBINS_DIR/clean_final_alignment.fasta}"
OUTDIR="${SEN_IQTREE_OUT:-$SEN_ROOT/iqtree_final}"
PREFIX="${SEN_IQTREE_PREFIX:-SSLAB_FINAL}"
CORES="${SEN_IQTREE_THREADS:-64}"
MEMORY="${SEN_IQTREE_MEMORY:-180G}"
MODEL="${SEN_IQTREE_MODEL:-MFP+ASC}"
BOOTSTRAPS="${SEN_IQTREE_BOOTSTRAPS:-1000}"

require_file "$INPUT_FILE"
mkdir -p "$OUTDIR"

echo "[INFO] IQ-TREE input: $INPUT_FILE"
echo "[INFO] Model: $MODEL | UFBoot: $BOOTSTRAPS | BNNI enabled"

iqtree2 -s "$INPUT_FILE"   -pre "$OUTDIR/$PREFIX"   -st DNA   -m "$MODEL"   -bb "$BOOTSTRAPS"   -bnni   -nt "$CORES"   -mem "$MEMORY"

echo "[INFO] Analysis complete: $OUTDIR/$PREFIX.treefile"
