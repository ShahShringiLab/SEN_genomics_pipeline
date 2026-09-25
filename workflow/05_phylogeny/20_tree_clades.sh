#!/bin/bash
set -euo pipefail

BASE_DIR="/home/samuelajulo/SENBio/Final"

TREE_DIR="${BASE_DIR}/iqtree_final"
TREEFILE="${TREE_DIR}/SSLAB_FINAL.treefile"

ITOL_DIR="${BASE_DIR}/ITOL"
CLADE1="${ITOL_DIR}/Clade1.txt"
CLADE2="${ITOL_DIR}/Clade2.txt"
CLADE3="${ITOL_DIR}/Clade3.txt"
UNASSIGNED="${ITOL_DIR}/Unassigned.txt"

tmp_all="$(mktemp)"
tmp_clades="$(mktemp)"
trap 'rm -f "$tmp_all" "$tmp_clades"' EXIT

echo "[INFO] Using treefile: $TREEFILE"
[[ -s "$TREEFILE" ]] || { echo "[ERROR] Missing/empty treefile: $TREEFILE" >&2; exit 1; }

# ----------------------------
# 1) Extract ALL SRR leaf IDs from the treefile (ONLY SRRxxxx)
# ----------------------------
grep -oE 'SRR[0-9]+' "$TREEFILE" | sort -u > "$tmp_all"

ALL_N=$(wc -l < "$tmp_all" | tr -d ' ')
echo "[INFO] SRRs found in tree: $ALL_N"

# ----------------------------
# 2) Build union of clade SRRs (sanitize: keep only SRR IDs)
# ----------------------------
for f in "$CLADE1" "$CLADE2" "$CLADE3"; do
  [[ -s "$f" ]] || { echo "[ERROR] Missing/empty clade file: $f" >&2; exit 1; }
done

cat "$CLADE1" "$CLADE2" "$CLADE3" \
  | tr -d '\r' \
  | grep -oE 'SRR[0-9]+' \
  | sort -u > "$tmp_clades"

CLADES_N=$(wc -l < "$tmp_clades" | tr -d ' ')
echo "[INFO] Unique SRRs across clades: $CLADES_N"

# ----------------------------
# 3) Unassigned = all_in_tree - clades_union
# ----------------------------
comm -23 "$tmp_all" "$tmp_clades" > "$UNASSIGNED"

UN_N=$(wc -l < "$UNASSIGNED" | tr -d ' ')
echo "[INFO] Unassigned SRRs written: $UN_N -> $UNASSIGNED"

# Optional: quick sanity sample
echo "[INFO] First 10 unassigned:"
head -n 10 "$UNASSIGNED" || true

