#!/usr/bin/env bash
# Replot only the three-panel Part A Hill diversity/composition figure.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"
export SEN_PARTA_SNAPSHOTS="${SEN_PARTA_SNAPSHOTS:-$ROOT}"
CONDA_BIN="$(command -v conda 2>/dev/null || true)"
if [[ -z "$CONDA_BIN" ]]; then
  for path in "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda"; do
    [[ ! -x "$path" ]] || { CONDA_BIN="$path"; break; }
  done
fi
[[ -n "$CONDA_BIN" ]] || { echo "[ERROR] Conda not found" >&2; exit 1; }
if [[ ! -s workflow/A2_Genotype/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd ]]; then
  python3 scripts/prepare_original_part_a2.py
fi
python3 scripts/build_part_a_diversity_figure.py
OUT="$ROOT/results/part_a/mdpi_label_corrections"
"$CONDA_BIN" run --no-capture-output -n sen_stats \
  Rscript "$OUT/run_diversity_figure_only.R" 2>&1 | tee "$OUT/diversity_figures.log"
rm -f "$OUT/QR_diversity_and_composition_side_by_side.png"
test -s "$OUT/QR_diversity_and_composition_side_by_side.tif"
echo "[PASS] TIFF only: $OUT/QR_diversity_and_composition_side_by_side.tif"
