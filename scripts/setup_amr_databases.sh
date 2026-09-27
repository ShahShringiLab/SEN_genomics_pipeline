#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEN_ROOT="${SEN_ROOT:-$REPO_ROOT}"
DB_ROOT="${SEN_AMR_DB_ROOT:-$SEN_ROOT/databases}"
AMR_ROOT="${SEN_AMRFINDER_DB_ROOT:-$DB_ROOT/amrfinderplus}"
PROV_DIR="${SEN_AMR_PROVENANCE_DIR:-$DB_ROOT/provenance}"
MARKER="$AMR_ROOT/.sen_verified"

mkdir -p "$AMR_ROOT" "$PROV_DIR"

command -v amrfinder >/dev/null 2>&1 || {
  echo "[ERROR] amrfinder not found." >&2
  exit 1
}
command -v amrfinder_update >/dev/null 2>&1 || {
  echo "[ERROR] amrfinder_update not found." >&2
  exit 1
}
command -v abricate >/dev/null 2>&1 || {
  echo "[ERROR] abricate not found." >&2
  exit 1
}

validate_amr_db() {
  [[ -e "$AMR_ROOT/latest" ]] || return 1
  local resolved
  resolved="$(readlink -f "$AMR_ROOT/latest" 2>/dev/null || true)"
  [[ -n "$resolved" && -d "$resolved" ]] || return 1
  [[ -s "$resolved/version.txt" ]] || return 1
  amrfinder -V --database "$resolved" >/dev/null 2>&1
}

if [[ -s "$MARKER" ]] && validate_amr_db; then
  echo "[SKIP] Frozen AMRFinderPlus database already validated."
else
  echo "[INFO] Bootstrapping AMRFinderPlus database into $AMR_ROOT ..."
  echo "[INFO] This update occurs only when the local frozen database is absent."
  amrfinder_update --database "$AMR_ROOT"

  if ! validate_amr_db; then
    rm -f "$MARKER"
    echo "[ERROR] AMRFinderPlus database bootstrap failed validation." >&2
    exit 1
  fi

  resolved="$(readlink -f "$AMR_ROOT/latest")"
  {
    echo "bootstrap_policy=freeze-on-first-success"
    echo "database_root=$AMR_ROOT"
    echo "resolved_database=$resolved"
    echo "resolved_database_name=$(basename "$resolved")"
    echo "bootstrap_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    amrfinder -V --database "$resolved" 2>&1 | sed 's/^/amrfinder_version_info=/'
  } > "$PROV_DIR/amrfinderplus.tsv"

  cp "$PROV_DIR/amrfinderplus.tsv" "$MARKER"
fi

AMR_DB="$(readlink -f "$AMR_ROOT/latest")"
echo "[INFO] Frozen AMRFinderPlus DB: $AMR_DB"
amrfinder -V --database "$AMR_DB" 2>&1 || true

echo "[INFO] Capturing bundled ABRicate database inventory..."
abricate --list > "$PROV_DIR/abricate_databases.tsv"

for db in resfinder vfdb plasmidfinder; do
  awk -v db="$db" 'BEGIN{FS="\t"} NR>1 && $1==db && $2+0>0 {ok=1} END{exit !ok}'     "$PROV_DIR/abricate_databases.tsv" || {
      echo "[ERROR] Required ABRicate database missing or empty: $db" >&2
      cat "$PROV_DIR/abricate_databases.tsv" >&2
      exit 1
    }
done

{
  echo "abricate_version=$(abricate --version 2>&1 | head -n1)"
  echo "database_policy=bundled-with-pinned-abricate-package-no-runtime-update"
  echo "inventory=$PROV_DIR/abricate_databases.tsv"
  echo "captured_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$PROV_DIR/abricate_provenance.txt"

echo "[PASS] AMR/ABRicate databases are ready."
echo "SEN_AMRFINDER_DB=$AMR_DB"
