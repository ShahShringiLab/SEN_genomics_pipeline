#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/conda_tool_inventory.tsv}"

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda is not available in PATH." >&2
  exit 1
}

TOOLS=(
  prefetch fasterq-dump fastp fastqc multiqc kraken2 seqkit
  SeqSero2_package.py skesa sistr mlst snippy snippy-core snpEff
  run_gubbins.py veryfasttree iqtree2 shovill spades.py prokka
  panaroo amrfinder abricate parallel python3 R
)

printf "environment\tenv_path\ttool\texecutable\n" > "$OUT"

mapfile -t ENVS < <(
  conda env list --json | python3 -c '
import json,sys
for p in json.load(sys.stdin).get("envs", []):
    print(p)
'
)

for env_path in "${ENVS[@]}"; do
  env_name="$(basename "$env_path")"
  echo "[INFO] Scanning $env_name"

  for tool in "${TOOLS[@]}"; do
    exe="$(
      conda run -p "$env_path" --no-capture-output         bash -lc "command -v '$tool' 2>/dev/null || true"         2>/dev/null | sed '/^[[:space:]]*$/d' | tail -n 1
    )"

    if [[ -n "$exe" ]]; then
      printf "%s\t%s\t%s\t%s\n"         "$env_name" "$env_path" "$tool" "$exe" >> "$OUT"
    fi
  done
done

echo "Wrote: $OUT"
echo
echo "Tools discovered per environment:"
awk -F'\t' 'NR>1 {n[$1]++} END {for (e in n) print e, n[e]}' "$OUT" | sort
