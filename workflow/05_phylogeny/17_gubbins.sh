#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd run_gubbins.py

ALN="${SEN_CORE_ALIGNMENT:-$SEN_ROOT/Snippy_output/senbio_core.full.aln}"
OUT_DIR="${SEN_GUBBINS_OUT:-$SEN_ROOT/gubbin}"
TMP_DIR="$OUT_DIR/tmp_gubbins"
THREADS="${SEN_GUBBINS_THREADS:-64}"
PREFIX="${SEN_GUBBINS_PREFIX:-senbio_res}"
ITERATIONS="${SEN_GUBBINS_ITERATIONS:-3}"

require_file "$ALN"
mkdir -p "$OUT_DIR" "$TMP_DIR"

if ! command -v veryfasttree >/dev/null 2>&1; then
    if [[ -n "${SEN_VERYFASTTREE_DIR:-}" && -x "$SEN_VERYFASTTREE_DIR/VeryFastTree" ]]; then
        export PATH="$SEN_VERYFASTTREE_DIR:$PATH"
    elif [[ -x "$SEN_ROOT/veryfasttree/build/VeryFastTree" ]]; then
        export PATH="$SEN_ROOT/veryfasttree/build:$PATH"
    else
        echo "[ERROR] veryfasttree not found on PATH." >&2
        exit 1
    fi
fi

export TMPDIR="$TMP_DIR"
cd "$OUT_DIR"

FILTERED="$OUT_DIR/${PREFIX}.filtered_polymorphic_sites.fasta"
if [[ -s "$FILTERED" ]]; then
  echo "[SKIP] Gubbins output already present: $FILTERED"
  exit 0
fi

echo "=== GUBBINS START: $(date) ==="
echo "Threads: $THREADS | Mode: --mar | Builder: veryfasttree"

run_gubbins.py   --threads "$THREADS"   --prefix "$PREFIX"   --tree-builder veryfasttree   --mar   --iterations "$ITERATIONS"   --verbose   "$ALN" 2>&1 | tee "$OUT_DIR/gubbins_run.log"

[[ -s "$FILTERED" ]] || {
  echo "[ERROR] Gubbins completed without expected filtered alignment: $FILTERED" >&2
  exit 1
}

echo "=== GUBBINS COMPLETE: $(date) ==="
