#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEN_ROOT="${SEN_ROOT:-$REPO_ROOT}"

# shellcheck disable=SC1091
source "$REPO_ROOT/config/database_sources.env"

DB_ROOT="${SEN_AMR_DB_ROOT:-$SEN_ROOT/databases}"
AMR_ROOT="${SEN_AMRFINDER_DB_ROOT:-$DB_ROOT/amrfinderplus}"
PROV_DIR="${SEN_AMR_PROVENANCE_DIR:-$DB_ROOT/provenance}"

AMR_MINOR="${SEN_AMRFINDER_DB_SOFTWARE_MINOR:?missing SEN_AMRFINDER_DB_SOFTWARE_MINOR}"
AMR_VERSION="${SEN_AMRFINDER_DB_VERSION:?missing SEN_AMRFINDER_DB_VERSION}"
AMR_BASE_URL="${SEN_AMRFINDER_DB_BASE_URL:?missing SEN_AMRFINDER_DB_BASE_URL}"
AMR_URL="$AMR_BASE_URL/$AMR_MINOR/$AMR_VERSION"
AMR_DIR="$AMR_ROOT/$AMR_VERSION"
LATEST_LINK="$AMR_ROOT/latest"
MARKER="$AMR_DIR/.sen_verified"

mkdir -p "$AMR_ROOT" "$PROV_DIR"

command -v amrfinder >/dev/null 2>&1 || {
  echo "[ERROR] amrfinder not found." >&2
  exit 1
}
command -v amrfinder_index >/dev/null 2>&1 || {
  echo "[ERROR] amrfinder_index not found." >&2
  exit 1
}
command -v abricate >/dev/null 2>&1 || {
  echo "[ERROR] abricate not found." >&2
  exit 1
}

fetch() {
  local url="$1" out="$2"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 5 --retry-delay 2 "$url" -o "$out"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
  else
    echo "[ERROR] curl or wget is required." >&2
    exit 1
  fi
}

validate_amr_db() {
  [[ -d "$AMR_DIR" && -s "$AMR_DIR/version.txt" ]] || return 1
  grep -q "$AMR_VERSION" "$AMR_DIR/version.txt" || return 1
  [[ -s "$AMR_DIR/AMRProt.fa.phr" ]] || return 1
  amrfinder -V --database "$AMR_DIR" 2>/dev/null | grep -q "Database version: $AMR_VERSION"
}

if [[ -s "$MARKER" ]] && validate_amr_db; then
  echo "[SKIP] Exact AMRFinderPlus database already validated: $AMR_VERSION"
else
  echo "[INFO] Bootstrapping exact AMRFinderPlus database release: $AMR_VERSION"
  echo "[INFO] Source: $AMR_URL/"

  rm -rf "$AMR_DIR"
  mkdir -p "$AMR_DIR"

  core_files=(
    AMR.LIB
    AMRProt.fa
    AMRProt-mutation.tsv
    AMRProt-suppress.tsv
    AMRProt-susceptible.fa
    AMRProt-susceptible.tsv
    AMR_CDS.fa
    database_format_version.txt
    fam.tsv
    taxgroup.tsv
    version.txt
    changes.txt
  )

  for name in "${core_files[@]}"; do
    echo "[INFO] Downloading $name"
    fetch "$AMR_URL/$name" "$AMR_DIR/$name"
  done

  while IFS= read -r taxgroup; do
    [[ -n "$taxgroup" ]] || continue
    for ext in fa tsv; do
      name="AMR_DNA-${taxgroup}.${ext}"
      echo "[INFO] Downloading $name"
      fetch "$AMR_URL/$name" "$AMR_DIR/$name"
    done
  done < <(
    awk '!/^#/ && NF>=3 && $3+0>0 {print $1}' "$AMR_DIR/taxgroup.tsv"
  )

  echo "[INFO] Indexing frozen AMRFinderPlus database..."
  amrfinder_index "$AMR_DIR"

  if ! validate_amr_db; then
    rm -f "$MARKER"
    echo "[ERROR] Exact AMRFinderPlus database failed validation." >&2
    exit 1
  fi

  printf 'database_version=%s\nsource=%s/\nverified_utc=%s\n'     "$AMR_VERSION" "$AMR_URL" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER"
fi

ln -sfn "$AMR_VERSION" "$LATEST_LINK"

{
  echo "database_policy=exact-release"
  echo "database_root=$AMR_ROOT"
  echo "resolved_database=$AMR_DIR"
  echo "database_version=$AMR_VERSION"
  echo "source_url=$AMR_URL/"
  echo "captured_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  amrfinder -V --database "$AMR_DIR" 2>&1 | sed 's/^/amrfinder_version_info=/'
} > "$PROV_DIR/amrfinderplus.tsv"

echo "[INFO] Capturing ABRicate database inventory..."
abricate --list > "$PROV_DIR/abricate_databases.tsv"

validate_abricate_row() {
  local db="$1" seqs="$2" date="$3"
  awk -F'\t' -v db="$db" -v seqs="$seqs" -v date="$date" '
    NR>1 && $1==db && $2==seqs && $4==date {ok=1}
    END{exit !ok}
  ' "$PROV_DIR/abricate_databases.tsv"
}

validate_abricate_row resfinder   "$SEN_ABRICATE_RESFINDER_SEQUENCES" "$SEN_ABRICATE_RESFINDER_DATE" || {
    echo "[ERROR] ResFinder inventory does not match pinned expectation." >&2
    exit 1
  }
validate_abricate_row vfdb   "$SEN_ABRICATE_VFDB_SEQUENCES" "$SEN_ABRICATE_VFDB_DATE" || {
    echo "[ERROR] VFDB inventory does not match pinned expectation." >&2
    exit 1
  }
validate_abricate_row plasmidfinder   "$SEN_ABRICATE_PLASMIDFINDER_SEQUENCES" "$SEN_ABRICATE_PLASMIDFINDER_DATE" || {
    echo "[ERROR] PlasmidFinder inventory does not match pinned expectation." >&2
    exit 1
  }

{
  echo "abricate_version=$(abricate --version 2>&1 | head -n1)"
  echo "database_policy=exact-bundled-inventory"
  echo "resfinder_sequences=$SEN_ABRICATE_RESFINDER_SEQUENCES"
  echo "resfinder_date=$SEN_ABRICATE_RESFINDER_DATE"
  echo "vfdb_sequences=$SEN_ABRICATE_VFDB_SEQUENCES"
  echo "vfdb_date=$SEN_ABRICATE_VFDB_DATE"
  echo "plasmidfinder_sequences=$SEN_ABRICATE_PLASMIDFINDER_SEQUENCES"
  echo "plasmidfinder_date=$SEN_ABRICATE_PLASMIDFINDER_DATE"
  echo "captured_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$PROV_DIR/abricate_provenance.txt"

echo "[PASS] Exact AMR/ABRicate database contracts validated."
echo "SEN_AMRFINDER_DB=$AMR_DIR"
