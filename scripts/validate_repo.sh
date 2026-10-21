#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

fail=0

echo "== SEN repository validation =="

echo
echo "[1/8] Hard-coded workstation paths"
if grep -RInE '/home/samuelajulo|SENBio/Final|~/miniforge3|conda activate'   workflow analysis config scripts   --include='*.sh' --include='*.py' --include='*.R'   --include='*.yaml' --include='*.yml'; then
  echo "[FAIL] Hard-coded workstation references remain."
  fail=1
else
  echo "[PASS] No hard-coded workstation references found."
fi

echo
echo "[2/8] Broken escaped shell-variable placeholders"
escaped_fail=0
while IFS= read -r f; do
  [[ "$f" == "scripts/validate_repo.sh" ]] && continue
  if grep -nF '\${' "$f"; then
    echo "[FAIL] literal escaped variable placeholder in $f"
    escaped_fail=1
  fi
done < <(find workflow scripts -type f -name '*.sh' | sort)

if [[ "$escaped_fail" -eq 0 ]]; then
  echo "[PASS] No broken escaped shell-variable placeholders found."
else
  fail=1
fi

echo
echo "[3/8] Bash syntax"
while IFS= read -r f; do
  if bash -n "$f"; then
    echo "[PASS] $f"
  else
    echo "[FAIL] $f"
    fail=1
  fi
done < <(find workflow scripts -type f -name '*.sh' | sort)

echo
echo "[4/8] Python syntax"
while IFS= read -r f; do
  if python3 -m py_compile "$f"; then
    echo "[PASS] $f"
  else
    echo "[FAIL] $f"
    fail=1
  fi
done < <(find workflow analysis scripts -type f -name '*.py' | sort)

echo
echo "[5/8] R syntax"
if command -v Rscript >/dev/null 2>&1; then
  R_CMD=(Rscript)
elif command -v conda >/dev/null 2>&1 && conda env list | awk '{print $1}' | grep -qx sen_stats; then
  R_CMD=(conda run --no-capture-output -n sen_stats Rscript)
else
  R_CMD=()
fi

if [[ "${#R_CMD[@]}" -eq 0 ]]; then
  echo "[SKIP] Rscript unavailable and sen_stats not installed; R parsing is checked during Part A runtime validation."
else
  while IFS= read -r f; do
    if "${R_CMD[@]}" -e 'parse(file=commandArgs(trailingOnly=TRUE)[1])' "$f" >/dev/null; then
      echo "[PASS] $f"
    else
      echo "[FAIL] $f"
      fail=1
    fi
  done < <(find analysis -type f -name '*.R' | sort)
fi

echo
echo "[6/8] Canonical metadata"
python3 - <<'PY'
import csv
from collections import Counter
from pathlib import Path

root = Path.cwd()
meta = root / "metadata" / "SEN_Genomes.csv"
clades = root / "metadata" / "final_clade_metadata.tsv"

with meta.open(newline="", encoding="utf-8-sig") as handle:
    m = list(csv.DictReader(handle))

with clades.open(newline="", encoding="utf-8-sig") as handle:
    c = list(csv.DictReader(handle, delimiter="\t"))

assert len(m) == 3434, f"SEN_Genomes.csv rows={len(m)}, expected 3434"
assert "Run" in m[0], "SEN_Genomes.csv missing Run column"
assert len({r["Run"] for r in m}) == 3434, "SEN_Genomes.csv Run values are not unique"

assert len(c) == 3307, f"final_clade_metadata rows={len(c)}, expected 3307"
assert {"Sample", "Lineage"}.issubset(c[0]), "clade metadata missing Sample/Lineage"
assert len({r["Sample"] for r in c}) == 3307, "final_clade_metadata Sample values are not unique"

counts = Counter(r["Lineage"] for r in c)
expected = Counter({
    "Clade 3": 1357,
    "Clade 2": 894,
    "Clade 1B": 759,
    "Clade 1A": 169,
    "Unassigned": 127,
    "Reference": 1,
})
assert counts == expected, f"Unexpected clade counts: {dict(counts)}"
print("[PASS] Canonical metadata dimensions and clade counts")
PY

echo
echo "[7/8] Runtime dependency installation guards"
runtime_install_fail=0
if grep -RInE 'install\.packages\(|BiocManager::install\(|pip[[:space:]]+install|conda[[:space:]]+install' \
  analysis workflow \
  --include='*.R' --include='*.py' --include='*.sh'; then
  echo "[FAIL] Runtime dependency installation found in analysis/workflow code."
  runtime_install_fail=1
else
  echo "[PASS] Analysis/workflow scripts do not install dependencies at runtime."
fi
if [[ "$runtime_install_fail" -ne 0 ]]; then
  fail=1
fi

echo
echo "[8/8] Git working tree"
git status --short

echo
if [[ "$fail" -eq 0 ]]; then
  echo "VALIDATION RESULT: PASS"
else
  echo "VALIDATION RESULT: FAIL"
  exit 1
fi
