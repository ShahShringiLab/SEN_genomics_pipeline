#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd mlst
require_cmd parallel
require_cmd nproc

IN_FASTA="${SEN_SISTR_ASSEMBLIES:-$SEN_ROOT/sistr_results_run/mini_assemblies}"
OUT_DIR="${SEN_MLST_OUT:-$SEN_ROOT/mlst_results_run}"
MASTER_REPORT="$OUT_DIR/mlst_master_report.csv"

mkdir -p "$OUT_DIR"

echo "=== MLST PIPELINE STARTED: $(date) ==="

CORES="${SEN_MLST_JOBS:-$(nproc)}"
echo "Utilizing $CORES parallel jobs for analysis..."

find "$IN_FASTA" -name "*.fasta" |   parallel -j "$CORES" -n 100 "mlst --threads 1 --csv {}"   > "$OUT_DIR/all_results.raw.csv"

echo "file,species,ST,aroC,dnaN,hemD,hisD,purE,sucA,thrA" > "$MASTER_REPORT"
cat "$OUT_DIR/all_results.raw.csv" >> "$MASTER_REPORT"
rm -f "$OUT_DIR/all_results.raw.csv"

echo "=== SUCCESS: Master MLST report created at $MASTER_REPORT ==="
echo "=== PIPELINE FINISHED: $(date) ==="
