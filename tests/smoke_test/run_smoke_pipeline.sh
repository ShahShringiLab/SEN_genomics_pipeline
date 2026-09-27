#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "=================================================="
echo " SEN reproducibility smoke pipeline"
echo "=================================================="
echo "[INFO] This runner bootstraps software/data dependencies as needed."
echo

python tests/smoke_test/validate_smoke.py
bash tests/smoke_test/run_reads_qc.sh
bash tests/smoke_test/run_kraken_typing.sh
bash tests/smoke_test/run_core_phylogeny.sh

echo
echo "=================================================="
echo " SEN SMOKE PIPELINE CURRENTLY IMPLEMENTED: PASS"
echo "=================================================="
echo "[INFO] Completed through Kraken2 + typing + core-SNP + Gubbins + IQ-TREE."
echo "[INFO] Pangenome and AMR stages remain to be appended after their"
echo "[INFO] database/bootstrap contracts are finalized."
