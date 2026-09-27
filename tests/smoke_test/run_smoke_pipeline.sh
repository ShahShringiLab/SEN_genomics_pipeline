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
bash tests/smoke_test/run_assembly_pangenome.sh

echo
echo "=================================================="
echo " SEN SMOKE PIPELINE CURRENTLY IMPLEMENTED: PASS"
echo "=================================================="
echo "[INFO] Completed through reads/QC + typing + phylogeny + assembly/pangenome."
echo "[INFO] AMR/virulence/plasmid screening remains to be appended after its"
echo "[INFO] database bootstrap contract is finalized."
