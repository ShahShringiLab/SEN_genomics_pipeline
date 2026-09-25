#!/bin/bash
set -e
set -o pipefail

# ========================================================
# GUBBINS - ULTIMATE STABLE PRODUCTION (3,400 TAXA)
# ========================================================

# 1. Safe Environment Activation
# We turn off -u temporarily to stop Conda from complaining about its own variables
set +u
source ~/miniforge3/etc/profile.d/conda.sh
conda activate gubbins_env
set -u # Turn strict mode back on for your variables

# 2. PATHS
WD="/home/samuelajulo/SENBio/Final"
ALN="$WD/Snippy_output/senbio_core.full.aln"
OUT_DIR="$WD/gubbin"
TMP_DIR="$OUT_DIR/tmp_gubbins"

mkdir -p "$OUT_DIR" "$TMP_DIR"
cd "$OUT_DIR"

# 3. VeryFastTree Check
# Ensuring the binary is found
if ! command -v veryfasttree &> /dev/null; then
    # If not in path, try looking in the common build location
    if [ -f "$WD/veryfasttree/build/VeryFastTree" ]; then
        export PATH="$WD/veryfasttree/build:$PATH"
    else
        echo "[ERROR] veryfasttree not found. Build it first!"
        exit 1
    fi
fi

# 4. ENVIRONMENT OVERRIDES
export TMPDIR="$TMP_DIR"

echo "=== GUBBINS START: $(date) ==="
echo "Threads: 64 | Mode: --mar | Builder: veryfasttree"

# 5. EXECUTION
# stdbuf ensures you see the progress in the log immediately
nohup stdbuf -oL -eL run_gubbins.py \
  --threads 64 \
  --prefix senbio_res \
  --tree-builder veryfasttree \
  --mar \
  --iterations 3 \
  --verbose \
  "$ALN" > gubbins_run.log 2>&1 &

echo $! > "$OUT_DIR/gubbins.pid"

echo "--------------------------------------------------------"
echo "PROCESS STARTED SUCCESSFULLY."
echo "PID: $(cat "$OUT_DIR/gubbins.pid")"
echo "Monitor: tail -f $OUT_DIR/gubbins_run.log"
echo "--------------------------------------------------------"
