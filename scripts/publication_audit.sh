#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "=================================================="
echo " SEN publication audit"
echo "=================================================="

bash scripts/validate_repo.sh
python scripts/validate_analysis_contracts.py

echo
echo "[AUDIT] Canonical final clade metadata"
python - <<'PY'
import pandas as pd
from pathlib import Path

p = Path("metadata/final_clade_metadata.tsv")
df = pd.read_csv(p, sep="\t")
expected = {
    "Clade 3": 1357,
    "Clade 2": 894,
    "Clade 1B": 759,
    "Clade 1A": 169,
    "Unassigned": 127,
    "Reference": 1,
}
counts = df["Lineage"].value_counts().to_dict()
assert len(df) == 3307, (len(df), 3307)
assert counts == expected, (counts, expected)
assert df["Sample"].nunique() == 3307
print("[PASS] 3307 unique tips and canonical clade counts match.")
PY

echo
echo "[AUDIT] Frozen AMR database contracts"
grep -q 'SEN_AMRFINDER_DB_VERSION="2026-08-07.1"' config/database_sources.env
grep -q 'SEN_ABRICATE_RESFINDER_SEQUENCES="3206"' config/database_sources.env
grep -q 'SEN_ABRICATE_VFDB_SEQUENCES="4592"' config/database_sources.env
grep -q 'SEN_ABRICATE_PLASMIDFINDER_SEQUENCES="488"' config/database_sources.env
echo "[PASS] AMR database versions/inventories are pinned."

echo
echo "=================================================="
echo " STATIC PUBLICATION AUDIT: PASS"
echo "=================================================="
echo "[NOTE] Full-study regenerated statistical outputs and unresolved historical"
echo "[NOTE] software provenance remain manual freeze gates; see docs/publication_audit.md."
