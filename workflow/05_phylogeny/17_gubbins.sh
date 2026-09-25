#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd run_gubbins.py
require_cmd nohup
require_cmd stdbuf

ALN="${SEN_CORE_ALIGNMENT:-$SEN_ROOT/Snippy_output/senbio_core.full.aln}"
OUT_DIR="${SEN_GUBBINS_OUT:-$SEN_ROOT/gubbin}"
TMP_DIR="$OUT_DIR/tmp_gubbins"
THREADS="${SEN_GUBBINS_THREADS:-64}"

require_file "$ALN"
mkdir -p "$OUT_DIR" "$TMP_DIR"

if ! command -v veryfasttree >/dev/null 2>&1; then
    if [[ -n "${SEN_VERYFASTTREE_DIR:-}" && -x "$SEN_VERYFASTTREE_DIR/VeryFastTree" ]]; then
        export PATH="$SEN_VERYFASTTREE_DIR:$PATH"
    elif [[ -x "$SEN_ROOT/veryfasttree/build/VeryFastTree" ]]; then
        export PATH="$SEN_ROOT/veryfasttree/build:$PATH"
    else
        echo "[ERROR] veryfasttree not found. Set SEN_VERYFASTTREE_DIR or put it on PATH." >&2
        exit 1
    fi
fi

export TMPDIR="$TMP_DIR"
cd "$OUT_DIR"

echo "=== GUBBINS START: $(date) ==="
echo "Threads: $THREADS | Mode: --mar | Builder: veryfasttree"

nohup stdbuf -oL -eL run_gubbins.py \
  --threads "$THREADS" \
  --prefix senbio_res \
  --tree-builder veryfasttree \
  --mar \
  --iterations 3 \
  --verbose \
  "$ALN" > gubbins_run.log 2>&1 &

echo $! > "$OUT_DIR/gubbins.pid"
echo "PID: $(cat "$OUT_DIR/gubbins.pid")"
echo "Monitor: tail -f $OUT_DIR/gubbins_run.log"
