#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd amrfinder
require_cmd parallel

ASSEMBLY_DIR="${SEN_AMR_ASSEMBLY_DIR:-${SEN_FINAL_CONTIGS_DIR:-$SEN_ROOT/Final_Contigs_Only}}"
OUT_DIR="${SEN_AMRFINDER_OUT:-$SEN_ROOT/AMRFinderplus}"
DB_DIR="${SEN_AMRFINDER_DB:-}"
JOBS="${SEN_AMRFINDER_JOBS:-4}"
THREADS="${SEN_AMRFINDER_THREADS:-4}"
ORG="${SEN_AMRFINDER_ORGANISM:-Salmonella}"

require_dir "$ASSEMBLY_DIR"
[[ -n "$DB_DIR" && -d "$DB_DIR" ]] || {
  echo "[ERROR] SEN_AMRFINDER_DB must point to a frozen AMRFinderPlus database." >&2
  exit 1
}

mkdir -p "$OUT_DIR/logs"
MANIFEST="$OUT_DIR/assemblies_manifest.tsv"
MASTER="$OUT_DIR/master_AMRFinder_report.tsv"
: > "$MANIFEST"

shopt -s nullglob
for fa in "$ASSEMBLY_DIR"/*.fasta "$ASSEMBLY_DIR"/*.fa "$ASSEMBLY_DIR"/*.fna; do
  [[ -s "$fa" ]] || continue
  sample="$(basename "$fa")"
  sample="${sample%%.*}"
  printf '%s\t%s\n' "$sample" "$fa" >> "$MANIFEST"
done

[[ -s "$MANIFEST" ]] || {
  echo "[ERROR] No assemblies found in $ASSEMBLY_DIR" >&2
  exit 1
}

run_one() {
  local sample="$1" fa="$2"
  local out="$OUT_DIR/${sample}.tsv"
  local log="$OUT_DIR/logs/${sample}.log"

  if [[ -s "$out" ]]; then
    echo "[SKIP] AMRFinderPlus $sample"
    return 0
  fi

  echo "[RUN] AMRFinderPlus $sample"
  amrfinder     -n "$fa"     --organism "$ORG"     --plus     --threads "$THREADS"     --database "$DB_DIR"     -o "$out" >"$log" 2>&1

  [[ -s "$out" ]] || {
    echo "[ERROR] AMRFinderPlus produced no report for $sample; see $log" >&2
    return 1
  }
  echo "[DONE] AMRFinderPlus $sample"
}
export -f run_one
export OUT_DIR THREADS ORG DB_DIR

parallel --colsep '\t' -j "$JOBS" --halt now,fail=1   run_one {1} {2} :::: "$MANIFEST"

first=""
while IFS=$'\t' read -r sample fa; do
  [[ -s "$OUT_DIR/${sample}.tsv" ]] || {
    echo "[ERROR] Missing AMRFinderPlus output for $sample" >&2
    exit 1
  }
  [[ -n "$first" ]] || first="$OUT_DIR/${sample}.tsv"
done < "$MANIFEST"

echo -e "Sample\t$(head -n1 "$first")" > "$MASTER"
while IFS=$'\t' read -r sample fa; do
  tail -n +2 "$OUT_DIR/${sample}.tsv" | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s,$0}' >> "$MASTER"
done < "$MANIFEST"

amrfinder -V --database "$DB_DIR" > "$OUT_DIR/database_version.txt" 2>&1

echo "[INFO] AMRFinderPlus complete: $MASTER"
