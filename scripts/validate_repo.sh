#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

fail=0

echo "== SEN repository validation =="

echo
echo "[1/5] Hard-coded workstation paths"
if grep -RInE '/home/samuelajulo|SENBio/Final|~/miniforge3|conda activate'   workflow analysis config   --include='*.sh' --include='*.py' --include='*.R'   --include='*.yaml' --include='*.yml'; then
  echo "[FAIL] Hard-coded workstation references remain."
  fail=1
else
  echo "[PASS] No hard-coded workstation references found."
fi

echo
echo "[2/5] Bash syntax"
while IFS= read -r f; do
  if bash -n "$f"; then
    echo "[PASS] $f"
  else
    echo "[FAIL] $f"
    fail=1
  fi
done < <(find workflow scripts -type f -name '*.sh' | sort)

echo
echo "[3/5] Python syntax"
while IFS= read -r f; do
  if python3 -m py_compile "$f"; then
    echo "[PASS] $f"
  else
    echo "[FAIL] $f"
    fail=1
  fi
done < <(find workflow analysis scripts -type f -name '*.py' | sort)

echo
echo "[4/5] Canonical metadata"
python3 - <<'PY'
from pathlib import Path
import pandas as pd

root = Path.cwd()
meta = root / "metadata" / "SEN_Genomes.csv"
clades = root / "metadata" / "final_clade_metadata.tsv"

m = pd.read_csv(meta, low_memory=False)
c = pd.read_csv(clades, sep=None, engine="python")

assert len(m) == 3434, f"SEN_Genomes.csv rows={len(m)}, expected 3434"
assert "Run" in m.columns, "SEN_Genomes.csv missing Run column"
assert len(c) == 3307, f"final_clade_metadata rows={len(c)}, expected 3307"
assert {"Sample", "Lineage"}.issubset(c.columns), "clade metadata missing Sample/Lineage"

counts = c["Lineage"].value_counts().to_dict()
expected = {
    "Clade 3": 1357,
    "Clade 2": 894,
    "Clade 1B": 759,
    "Clade 1A": 169,
    "Unassigned": 127,
    "Reference": 1,
}
assert counts == expected, f"Unexpected clade counts: {counts}"
print("[PASS] Canonical metadata dimensions and clade counts")
PY

echo
echo "[5/5] Git working tree"
git status --short

echo
if [[ "$fail" -eq 0 ]]; then
  echo "VALIDATION RESULT: PASS"
else
  echo "VALIDATION RESULT: FAIL"
  exit 1
fi
