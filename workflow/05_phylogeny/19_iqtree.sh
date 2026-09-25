#!/bin/bash
set -euo pipefail

# ==============================================================================
# IQ-TREE 2: CLEANED FINAL TREE (3,307 TAXA)
# ==============================================================================

INPUT_FILE="gubbin/clean_final_alignment.fasta"
OUTDIR="iqtree"
mkdir -p "$OUTDIR"

CORES=64
MEMORY="180G"

echo "Starting IQ-TREE 2: Model Selection + ML Tree + 1000 UFBoot..."

# Using MFP+ASC:
# 1. MFP = ModelFinder Plus (picks best model automatically)
# 2. +ASC = Ascertainment Bias Correction (REQUIRED for SNP alignments)
# 3. -bb 1000 = Ultrafast Bootstrap
# 4. -bnni = Optimizes tree topology to be more accurate

iqtree2 -s "$INPUT_FILE" \
  -pre "${OUTDIR}/SSLAB_FINAL" \
  -st DNA \
  -m MFP+ASC \
  -bb 1000 \
  -bnni \
  -nt "$CORES" \
  -mem "$MEMORY"

echo "Analysis complete. Check ${OUTDIR}/SSLAB_FINAL.treefile"
