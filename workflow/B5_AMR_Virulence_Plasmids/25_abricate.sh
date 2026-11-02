#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd abricate
require_cmd parallel

ASSEMBLY_DIR="${SEN_ABRICATE_ASSEMBLY_DIR:-${SEN_FINAL_CONTIGS_DIR:-$SEN_ROOT/Final_Contigs_Only}}"
OUT_ROOT="${SEN_ABRICATE_OUT:-$SEN_ROOT/abricate_results}"
JOBS="${SEN_ABRICATE_JOBS:-4}"
THREADS="${SEN_ABRICATE_THREADS:-2}"
MINCOV="${SEN_ABRICATE_MINCOV:-80}"
MINID="${SEN_ABRICATE_MINID:-90}"
DATABASES=(resfinder vfdb plasmidfinder)

require_dir "$ASSEMBLY_DIR"
mkdir -p "$OUT_ROOT"

MANIFEST="$OUT_ROOT/assemblies_manifest.tsv"
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

abricate --list > "$OUT_ROOT/database_inventory.tsv"

for DB in "${DATABASES[@]}"; do
  if ! awk -v db="$DB" 'BEGIN{FS="\t"} NR>1 && $1==db && $2+0>0 {ok=1} END{exit !ok}'       "$OUT_ROOT/database_inventory.tsv"; then
    echo "[ERROR] Required ABRicate DB unavailable: $DB" >&2
    exit 1
  fi

  DB_OUT="$OUT_ROOT/$DB"
  mkdir -p "$DB_OUT"

  run_one() {
    local sample="$1" fa="$2"
    local out="$DB_OUT/${sample}_${DB}.tsv"
    if [[ -s "$out" ]]; then
      echo "[SKIP] ABRicate $DB $sample"
      return 0
    fi
    echo "[RUN] ABRicate $DB $sample"
    abricate       --db "$DB"       --mincov "$MINCOV"       --minid "$MINID"       --threads "$THREADS"       "$fa" > "$out"
    [[ -s "$out" ]] || {
      echo "[ERROR] ABRicate produced no report for $DB/$sample" >&2
      return 1
    }
  }
  export -f run_one
  export DB_OUT DB MINCOV MINID THREADS

  parallel --colsep '\t' -j "$JOBS" --halt now,fail=1     run_one {1} {2} :::: "$MANIFEST"

  FIRST="$(find "$DB_OUT" -maxdepth 1 -type f -name "*_${DB}.tsv" -size +0c | sort | head -n1)"
  [[ -n "$FIRST" ]] || {
    echo "[ERROR] No ABRicate reports found for $DB" >&2
    exit 1
  }

  MASTER="$DB_OUT/master_${DB}_report.tsv"
  echo -e "Sample\t$(head -n1 "$FIRST")" > "$MASTER"

  while IFS=$'\t' read -r sample fa; do
    report="$DB_OUT/${sample}_${DB}.tsv"
    [[ -s "$report" ]] || {
      echo "[ERROR] Missing ABRicate report: $report" >&2
      exit 1
    }
    tail -n +2 "$report" | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s,$0}' >> "$MASTER"
  done < "$MANIFEST"

  echo "[INFO] ABRicate $DB master: $MASTER"
done

echo "[INFO] ABRicate complete with mincov=$MINCOV minid=$MINID"
