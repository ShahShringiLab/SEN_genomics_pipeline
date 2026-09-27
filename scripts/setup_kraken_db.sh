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
: "${SEN_KRAKEN_DB_MD5_URL:=https://genome-idx.s3.amazonaws.com/kraken/standard_20260626/standard.md5}"

SEN_ROOT="${SEN_ROOT:-$REPO_ROOT}"
DB_PARENT="${SEN_DATABASE_DIR:-$SEN_ROOT/databases}"
KRAKEN_PARENT="$DB_PARENT/kraken2"
DB_DIR="${SEN_KRAKEN_DB:-$KRAKEN_PARENT/k2_standard_${SEN_KRAKEN_DB_DATE}}"
CACHE_DIR="${SEN_DB_CACHE_DIR:-$DB_PARENT/downloads}"
ARCHIVE="$CACHE_DIR/$SEN_KRAKEN_DB_ARCHIVE"
MD5_FILE="$CACHE_DIR/standard_${SEN_KRAKEN_DB_DATE}.md5"
VERIFY_MARKER="$DB_DIR/.sen_verified"
EXTRACT_MARKER="$DB_DIR/.sen_extracted"

mkdir -p "$KRAKEN_PARENT" "$CACHE_DIR" "$DB_DIR"

core_db_present() {
  [[ -s "$DB_DIR/hash.k2d" && -s "$DB_DIR/opts.k2d" && -s "$DB_DIR/taxo.k2d" ]]
}

is_verified_db() {
  core_db_present && [[ -s "$VERIFY_MARKER" ]]
}

if is_verified_db; then
  echo "[INFO] Verified Kraken2 database already present: $DB_DIR"
  exit 0
fi

if command -v aria2c >/dev/null 2>&1; then
  DOWNLOAD_TOOL="aria2c"
elif command -v wget >/dev/null 2>&1; then
  DOWNLOAD_TOOL="wget"
else
  echo "[ERROR] aria2c or wget is required for database download." >&2
  exit 1
fi
command -v tar >/dev/null 2>&1 || {
  echo "[ERROR] tar is required to unpack the Kraken2 database." >&2
  exit 1
}
if command -v pigz >/dev/null 2>&1; then
  EXTRACT_TOOL="pigz"
else
  EXTRACT_TOOL="gzip"
fi
command -v md5sum >/dev/null 2>&1 || {
  echo "[ERROR] md5sum is required to verify the Kraken2 database." >&2
  exit 1
}

echo "=================================================="
echo " Kraken2 database bootstrap"
echo "=================================================="
echo "[INFO] Snapshot:   Standard $SEN_KRAKEN_DB_DATE"
echo "[INFO] URL:        $SEN_KRAKEN_DB_URL"
echo "[INFO] Install to: $DB_DIR"
echo "[INFO] Cache:      $CACHE_DIR"
echo "[INFO] Downloader: $DOWNLOAD_TOOL"
echo

if [[ ! -s "$MD5_FILE" ]]; then
  echo "[INFO] Downloading official checksum manifest..."
  if [[ "$DOWNLOAD_TOOL" == "aria2c" ]]; then
    aria2c -c -x 4 -s 4 --file-allocation=none       -d "$CACHE_DIR" -o "$(basename "$MD5_FILE")" "$SEN_KRAKEN_DB_MD5_URL"
  else
    wget -O "$MD5_FILE" "$SEN_KRAKEN_DB_MD5_URL"
  fi
fi

if [[ ! -s "$ARCHIVE" ]]; then
  # Full Standard 2026-06-26 is ~79.6 GB compressed and ~103 GB unpacked.
  # Keep a safety margin because download and extraction coexist temporarily.
  avail_kb="$(df -Pk "$DB_PARENT" | awk 'NR==2 {print $4}')"
  required_kb=$((190 * 1024 * 1024))
  if [[ "$avail_kb" -lt "$required_kb" ]]; then
    echo "[ERROR] Less than ~190 GB free under $DB_PARENT." >&2
    echo "[ERROR] The pinned full Kraken2 Standard database needs substantial temporary space." >&2
    exit 1
  fi

  echo "[INFO] Downloading Kraken2 Standard archive (~80 GB compressed)..."
  if [[ "$DOWNLOAD_TOOL" == "aria2c" ]]; then
    ARIA2_CONNECTIONS="${SEN_ARIA2_CONNECTIONS:-16}"
    aria2c -c       -x "$ARIA2_CONNECTIONS"       -s "$ARIA2_CONNECTIONS"       -k 4M       --file-allocation=none       --summary-interval=10       -d "$CACHE_DIR"       -o "$SEN_KRAKEN_DB_ARCHIVE"       "$SEN_KRAKEN_DB_URL"
  else
    echo "[WARN] aria2c unavailable; falling back to single-connection wget."
    wget -c -O "$ARCHIVE" "$SEN_KRAKEN_DB_URL"
  fi
else
  echo "[SKIP] Archive already present: $ARCHIVE"
fi

echo "[INFO] Verifying downloaded archive against official manifest..."
archive_line="$(awk -v a="$SEN_KRAKEN_DB_ARCHIVE" '$2==a || $2=="*"a {print; exit}' "$MD5_FILE")"
if [[ -z "$archive_line" ]]; then
  echo "[ERROR] Official checksum manifest does not contain $SEN_KRAKEN_DB_ARCHIVE." >&2
  exit 1
fi
(
  cd "$CACHE_DIR"
  printf '%s\n' "$archive_line" | md5sum -c -
)

if [[ -s "$EXTRACT_MARKER" ]]; then
  echo "[SKIP] Extraction marker present; reusing extracted Kraken2 database."
else
  echo "[INFO] Extracting database..."
  if [[ "$EXTRACT_TOOL" == "pigz" ]]; then
    EXTRACT_THREADS="${SEN_PIGZ_EXTRACT_THREADS:-16}"
    echo "[INFO] Parallel decompression: pigz -p $EXTRACT_THREADS"
    pigz -p "$EXTRACT_THREADS" -dc "$ARCHIVE" | tar -xf - -C "$DB_DIR"
  else
    echo "[WARN] pigz unavailable; falling back to single-threaded gzip extraction."
    tar -xzf "$ARCHIVE" -C "$DB_DIR"
  fi
  touch "$EXTRACT_MARKER"
fi

echo "[INFO] Verifying extracted Kraken2 database files..."
internal_manifest="$DB_DIR/.standard_internal.md5"
awk -v a="$SEN_KRAKEN_DB_ARCHIVE" '$2!=a && $2!="*"a {print}' "$MD5_FILE" > "$internal_manifest"
if ! (
  cd "$DB_DIR"
  md5sum -c "$(basename "$internal_manifest")"
); then
  rm -f "$EXTRACT_MARKER" "$VERIFY_MARKER" "$internal_manifest"
  echo "[ERROR] Extracted Kraken2 database verification failed." >&2
  exit 1
fi
rm -f "$internal_manifest"

if ! core_db_present; then
  echo "[ERROR] Required Kraken2 database files are missing after verification." >&2
  exit 1
fi

cat > "$DB_DIR/SEN_DATABASE_PROVENANCE.txt" <<EOF
database=Kraken2 Standard
snapshot_date=$SEN_KRAKEN_DB_DATE
archive=$SEN_KRAKEN_DB_ARCHIVE
source_url=$SEN_KRAKEN_DB_URL
checksum_url=$SEN_KRAKEN_DB_MD5_URL
verified_archive=md5
verified_extracted_files=md5
installed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

printf 'snapshot=%s\nverified_utc=%s\n'   "$SEN_KRAKEN_DB_DATE" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$VERIFY_MARKER"

if [[ "${SEN_KEEP_DB_ARCHIVE:-0}" != "1" ]]; then
  rm -f "$ARCHIVE"
fi

echo "[PASS] Kraken2 database ready and verified: $DB_DIR"
