#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] Conda not found. Install Miniforge or initialize conda first." >&2
  exit 1
fi

if ! conda env list | awk '{print $1}' | grep -Fxq sen_figure5b; then
  echo "[INFO] Creating dedicated R environment sen_figure5b ..."
  conda env create -f environments/figure5b.yaml
fi

echo "[INFO] Running Figure 5B plot (trimmed at the observed data range) ..."
conda run --no-capture-output -n sen_figure5b \
  Rscript scripts/figure5b_violin_zero_bounded.R "$@"

echo "[PASS] Outputs: $ROOT/results/figure5b_zero_bounded/"
