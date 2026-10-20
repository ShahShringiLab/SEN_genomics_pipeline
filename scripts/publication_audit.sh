#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "=================================================="
echo " SEN publication audit"
echo "=================================================="

bash scripts/validate_repo.sh
python3 scripts/validate_analysis_contracts.py

echo
echo "[AUDIT] Part A canonical inputs and Part A -> Part B handoff"
python3 scripts/validate_part_a_inputs.py

echo
echo "[AUDIT] Canonical final clade metadata"
python3 - <<'PY'
import csv
from collections import Counter

with open("metadata/final_clade_metadata.tsv", newline="", encoding="utf-8-sig") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

expected = Counter({
    "Clade 3": 1357,
    "Clade 2": 894,
    "Clade 1B": 759,
    "Clade 1A": 169,
    "Unassigned": 127,
    "Reference": 1,
})
counts = Counter(r["Lineage"] for r in rows)
assert len(rows) == 3307, (len(rows), 3307)
assert counts == expected, (dict(counts), dict(expected))
assert len({r["Sample"] for r in rows}) == 3307
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
echo "[NOTE] Runtime Part A validation: bash scripts/run_part_a.sh"
echo "[NOTE] Manuscript-scale Part B execution: docs/FULL_RUN_FRESH_LINUX.md"
