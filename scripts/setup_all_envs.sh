#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda not found. Install Miniforge/Conda first." >&2
  exit 1
}

create_or_verify() {
  local yaml="$1"
  local name
  name="$(awk '$1=="name:" {print $2; exit}' "$yaml")"
  [[ -n "$name" ]] || { echo "[ERROR] no env name in $yaml" >&2; exit 1; }

  if conda env list | awk '{print $1}' | grep -qx "$name"; then
    echo "[SKIP] $name already exists"
  else
    echo "[CREATE] $name <- $yaml"
    if conda env create --help 2>&1 | grep -q -- '--solver'; then
      conda env create --solver=libmamba -f "$yaml"
    else
      conda env create -f "$yaml"
    fi
  fi
}

for yaml in environments/*.yaml; do
  create_or_verify "$yaml"
done

echo "[PASS] All declared SEN environments are present."
