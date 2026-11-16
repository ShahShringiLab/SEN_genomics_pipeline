#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Miniforge can be installed without its conda entry point being on PATH in WSL.
CONDA_BIN="$(command -v conda 2>/dev/null || true)"
if [[ -z "$CONDA_BIN" ]]; then
  for candidate in \
    "$HOME/miniforge3/bin/conda" \
    "$HOME/miniforge3/condabin/conda" \
    "$HOME/miniconda3/bin/conda" \
    "$HOME/miniconda3/condabin/conda"; do
    if [[ -x "$candidate" ]]; then
      CONDA_BIN="$candidate"
      break
    fi
  done
fi

if [[ -z "$CONDA_BIN" ]]; then
  echo "[ERROR] Conda was not located. Check the Miniforge install directory." >&2
  exit 1
fi
echo "[INFO] Conda executable: $CONDA_BIN"

if ! "$CONDA_BIN" env list | awk '{print $1}' | grep -Fxq sen_figure5b; then
  echo "[INFO] Creating dedicated R environment sen_figure5b ..."
  "$CONDA_BIN" env create -f environments/figure5b.yaml
fi

echo "[INFO] Running Figure 5B plot (trimmed at observed data range) ..."
"$CONDA_BIN" run --no-capture-output -n sen_figure5b \
  Rscript scripts/figure5b_violin_zero_bounded.R "$@"

echo "[PASS] Outputs: $ROOT/results/figure5b_zero_bounded/"
