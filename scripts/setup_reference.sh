#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEN_ROOT="${SEN_ROOT:-$REPO_ROOT}"
REF_DIR="${SEN_REFERENCE_DIR:-$SEN_ROOT/reference}"
ACCESSION="${SEN_REFERENCE_ACCESSION:-NC_011294.1}"
REF_FNA="${SEN_REFERENCE_FNA:-$REF_DIR/reference.fna}"
REF_GBK="${SEN_REFERENCE_GBK:-$REF_DIR/reference.gbk}"
MARKER="$REF_DIR/.sen_reference_verified"

FASTA_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACCESSION}&rettype=fasta&retmode=text"
GBK_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACCESSION}&rettype=gbwithparts&retmode=text"

mkdir -p "$REF_DIR"

fetch() {
  local url="$1" out="$2"
  if command -v curl >/dev/null 2>&1; then
    curl -fL --retry 5 --retry-delay 2 "$url" -o "$out"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$out" "$url"
  else
    echo "[ERROR] curl or wget is required to download the reference." >&2
    exit 1
  fi
}

validate_reference() {
  [[ -s "$REF_FNA" && -s "$REF_GBK" ]] || return 1
  grep -q "$ACCESSION" "$REF_FNA" || return 1
  grep -q "$ACCESSION" "$REF_GBK" || return 1

  local len
  len="$(awk '!/^>/ {gsub(/[[:space:]]/,""); n+=length($0)} END{print n+0}' "$REF_FNA")"
  # P125109 chromosome is ~4.7 Mb; broad bounds catch HTML/error payloads and truncation.
  [[ "$len" -ge 4500000 && "$len" -le 5000000 ]] || return 1
}

if [[ -s "$MARKER" ]] && validate_reference; then
  echo "[SKIP] Verified reference already present: $ACCESSION"
  exit 0
fi

echo "[INFO] Bootstrapping reference $ACCESSION from NCBI..."
tmp_fna="$REF_FNA.tmp"
tmp_gbk="$REF_GBK.tmp"
trap 'rm -f "$tmp_fna" "$tmp_gbk"' EXIT

fetch "$FASTA_URL" "$tmp_fna"
fetch "$GBK_URL" "$tmp_gbk"
mv "$tmp_fna" "$REF_FNA"
mv "$tmp_gbk" "$REF_GBK"

if ! validate_reference; then
  rm -f "$MARKER"
  echo "[ERROR] Downloaded reference failed validation." >&2
  exit 1
fi

cat > "$REF_DIR/SEN_REFERENCE_PROVENANCE.txt" <<EOF
organism=Salmonella enterica subsp. enterica serovar Enteritidis str. P125109
accession=$ACCESSION
fasta_source=$FASTA_URL
genbank_source=$GBK_URL
downloaded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF

printf 'accession=%s\nverified_utc=%s\n'   "$ACCESSION" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$MARKER"

echo "[PASS] Reference ready: $REF_FNA"
