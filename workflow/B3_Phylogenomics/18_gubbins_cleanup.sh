#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

GUBBINS_DIR="${SEN_GUBBINS_OUT:-$SEN_ROOT/gubbin}"
INPUT_ALIGNMENT="${SEN_GUBBINS_FILTERED_ALIGNMENT:-$GUBBINS_DIR/senbio_res.filtered_polymorphic_sites.fasta}"
OUTPUT_ALIGNMENT="${SEN_CLEAN_ALIGNMENT:-$GUBBINS_DIR/clean_final_alignment.fasta}"
EXCLUSION_FILE="${SEN_FINAL_TREE_EXCLUSIONS:-$REPO_ROOT/config/final_tree_exclusions.txt}"

require_file "$INPUT_ALIGNMENT"
require_file "$EXCLUSION_FILE"

tmp_exclusions="$(mktemp)"
trap 'rm -f "$tmp_exclusions"' EXIT

grep -vE '^[[:space:]]*(#|$)' "$EXCLUSION_FILE" | tr -d '\r' | sort -u > "$tmp_exclusions"

awk -v ex="$tmp_exclusions" '
BEGIN {
  while ((getline line < ex) > 0) {
    excluded[">" line] = 1
  }
}
/^>/ {
  key=$1
  skip=(key in excluded)
}
!skip
' "$INPUT_ALIGNMENT" > "$OUTPUT_ALIGNMENT"

orig=$(grep -c '^>' "$INPUT_ALIGNMENT" || true)
clean=$(grep -c '^>' "$OUTPUT_ALIGNMENT" || true)
removed=$((orig-clean))

echo "[INFO] Original alignment sequences: $orig"
echo "[INFO] Cleaned alignment sequences:  $clean"
echo "[INFO] Removed:                      $removed"
echo "[INFO] Output: $OUTPUT_ALIGNMENT"

if [[ "$clean" -eq 0 ]]; then
  echo "[ERROR] Cleaned alignment is empty." >&2
  exit 1
fi
