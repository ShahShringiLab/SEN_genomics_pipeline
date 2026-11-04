#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

TREEFILE="${SEN_TREEFILE:-${SEN_IQTREE_OUT:-$SEN_ROOT/iqtree_final}/SSLAB_FINAL.treefile}"
CLADE_META="${SEN_CLADE_METADATA:-$SEN_ROOT/metadata/final_clade_metadata.tsv}"

require_file "$TREEFILE"
require_file "$CLADE_META"

tmp_tree="$(mktemp)"
tmp_meta="$(mktemp)"
tmp_dup="$(mktemp)"
trap 'rm -f "$tmp_tree" "$tmp_meta" "$tmp_dup"' EXIT

grep -oE 'SRR[0-9]+' "$TREEFILE" | sort -u > "$tmp_tree"

awk -F'\t' '
  NR==1 {next}
  $1 ~ /^SRR[0-9]+$/ {
    if ($2!="Clade 1A" && $2!="Clade 1B" && $2!="Clade 2" && $2!="Clade 3" && $2!="Unassigned") {
      print "[ERROR] invalid lineage for " $1 ": " $2 > "/dev/stderr"
      bad=1
    }
    print $1
  }
  END {if (bad) exit 2}
' "$CLADE_META" | sort > "$tmp_meta"

awk -F'\t' 'NR>1 && $1 ~ /^SRR[0-9]+$/ {n[$1]++} END {for (x in n) if (n[x]>1) print x}'   "$CLADE_META" | sort > "$tmp_dup"

if [[ -s "$tmp_dup" ]]; then
  echo "[ERROR] Duplicate SRR assignments in canonical clade metadata:" >&2
  cat "$tmp_dup" >&2
  exit 1
fi

tree_n="$(wc -l < "$tmp_tree" | tr -d ' ')"
meta_n="$(wc -l < "$tmp_meta" | tr -d ' ')"

echo "[INFO] Tree SRR tips:            $tree_n"
echo "[INFO] Canonical metadata SRRs:  $meta_n"

if ! diff -u "$tmp_tree" "$tmp_meta"; then
  echo "[ERROR] Canonical clade metadata SRR set does not match the IQ-TREE SRR set." >&2
  exit 1
fi

echo "[INFO] Canonical lineage counts:"
awk -F'\t' 'NR>1 {n[$2]++} END {for (x in n) print x "\t" n[x]}' "$CLADE_META" | sort

echo "[PASS] Final IQ-TREE SRR set exactly matches metadata/final_clade_metadata.tsv."
echo "[INFO] Clade assignments are treated as frozen study metadata; this script does not re-derive them from historical ITOL lists."
