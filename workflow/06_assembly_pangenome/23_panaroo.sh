#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd panaroo
require_cmd mafft

GFF_DIR="${SEN_PROKKA_OUT:-$SEN_ROOT/prokka_annotations}"
OUT_DIR="${SEN_PANAROO_OUT:-$SEN_ROOT/Panaroo_Run/results}"
THREADS="${SEN_PANAROO_THREADS:-$(nproc)}"
CLEAN_MODE="${SEN_PANAROO_CLEAN_MODE:-strict}"

require_dir "$GFF_DIR"
mkdir -p "$OUT_DIR"

mapfile -t GFFS < <(find "$GFF_DIR" -maxdepth 1 -type f -name "*.gff" | sort)
if [[ ${#GFFS[@]} -eq 0 ]]; then
  echo "[ERROR] No GFF files found in $GFF_DIR" >&2
  exit 1
fi

echo "[INFO] Running Panaroo on ${#GFFS[@]} GFF files"
echo "[INFO] clean-mode=$CLEAN_MODE threads=$THREADS"

panaroo   -i "${GFFS[@]}"   -o "$OUT_DIR"   --clean-mode "$CLEAN_MODE"   --remove-invalid-genes   -t "$THREADS"   -a core   --aligner mafft

echo "[INFO] Panaroo complete: $OUT_DIR"
