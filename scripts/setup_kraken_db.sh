#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SOURCE_FILE="$REPO_ROOT/config/database_sources.env"

[[ -f "$SOURCE_FILE" ]] || {
  echo "[ERROR] Missing database source manifest: $SOURCE_FILE" >&2
  exit 1
}

# shellcheck disable=SC1090
source "$SOURCE_FILE"

: "${SEN_KRAKEN_DB_DATE:=20260626}"
: "${SEN_KRAKEN_DB_ARCHIVE:=k2_standard_20260626.tar.gz}"
: "${SEN_KRAKEN_DB_URL:=https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260626.tar.gz}"
: "${SEN_KRAKEN_DB_MD5_URL:=https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260626.tar.gz.md5}"

SEN_ROOT="${SEN_ROOT:-$REPO_ROOT}"
DB_PARENT="${SEN_DATABASE_DIR:-$SEN_ROOT/databases}"
KRAKEN_PARENT="$DB_PARENT/kraken2"
DB_DIR="${SEN_KRAKEN_DB:-$KRAKEN_PARENT/k2_standard_${SEN_KRAKEN_DB_DATE}}"
CACHE_DIR="${SEN_DB_CACHE_DIR:-$DB_PARENT/downloads}"
ARCHIVE="$CACHE_DIR/$SEN_KRAKEN_DB_ARCHIVE"
MD5_FILE="$ARCHIVE.md5"

mkdir -p "$KRAKEN_PARENT" "$CACHE_DIR" "$DB_DIR"

is_valid_db() {
  [[ -s "$DB_DIR/hash.k2d" && -s "$DB_DIR/opts.k2d" && -s "$DB_DIR/taxo.k2d" ]]
}

if is_valid_db; then
  echo "[INFO] Kraken2 database already present: $DB_DIR"
  exit 0
fi

command -v wget >/dev/null 2>&1 || {
  echo "[ERROR] wget is required for resumable database download." >&2
  exit 1
}
command -v tar >/dev/null 2>&1 || {
  echo "[ERROR] tar is required to unpack the Kraken2 database." >&2
  exit 1
}
command -v md5sum >/dev/null 2>&1 || {
  echo "[ERROR] md5sum is required to verify the Kraken2 database archive." >&2
  exit 1
}

# Full Standard 2026-06-26 is ~79.6 GB compressed and ~103 GB unpacked.
# Keep a safety margin because download and extraction coexist temporarily.
avail_kb="$(df -Pk "$DB_PARENT" | awk 'NR==2 {print $4}')"
required_kb=$((190 * 1024 * 1024))
if [[ "$avail_kb" -lt "$required_kb" ]]; then
  echo "[ERROR] Less than ~190 GB free under $DB_PARENT." >&2
  echo "[ERROR] The pinned full Kraken2 Standard database needs substantial temporary space." >&2
  exit 1
fi

echo "=================================================="
echo " Kraken2 database bootstrap"
echo "=================================================="
echo "[INFO] Snapshot:   Standard $SEN_KRAKEN_DB_DATE"
echo "[INFO] URL:        $SEN_KRAKEN_DB_URL"
echo "[INFO] Install to: $DB_DIR"
echo "[INFO] Cache:      $CACHE_DIR"
echo
echo "[INFO] This is a large download (~80 GB compressed). wget -c will resume."
echo

wget -c -O "$ARCHIVE" "$SEN_KRAKEN_DB_URL"
wget -O "$MD5_FILE" "$SEN_KRAKEN_DB_MD5_URL"

echo "[INFO] Verifying archive checksum..."
(
  cd "$CACHE_DIR"
  md5sum -c "$(basename "$MD5_FILE")"
)

echo "[INFO] Extracting database..."
tar -xzf "$ARCHIVE" -C "$DB_DIR"

if ! is_valid_db; then
  echo "[ERROR] Extraction completed but required Kraken2 files are missing." >&2
  exit 1
fi

cat > "$DB_DIR/SEN_DATABASE_PROVENANCE.txt" <<EOF
database=Kraken2 Standard
snapshot_date=$SEN_KRAKEN_DB_DATE
archive=$SEN_KRAKEN_DB_ARCHIVE
source_url=$SEN_KRAKEN_DB_URL
checksum_url=$SEN_KRAKEN_DB_MD5_URL
installed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

if [[ "${SEN_KEEP_DB_ARCHIVE:-0}" != "1" ]]; then
  rm -f "$ARCHIVE"
fi

echo "[PASS] Kraken2 database ready: $DB_DIR"
