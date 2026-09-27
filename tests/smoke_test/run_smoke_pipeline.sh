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
bash tests/smoke_test/run_amr_vf_plasmid.sh

echo
echo "=================================================="
echo " SEN END-TO-END SMOKE PIPELINE: PASS"
echo "=================================================="
echo "[INFO] Completed reads/QC, typing, core-SNP phylogeny, assembly/pangenome,"
echo "[INFO] and AMR/virulence/plasmid screening."
