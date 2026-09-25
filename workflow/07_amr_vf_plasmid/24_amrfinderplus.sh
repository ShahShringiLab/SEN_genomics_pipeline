#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

# =========================================================
# 06_amrfinderplus.sh  (ALL_ASSEMBLIES version)
# - Archives current AMRFinderplus outputs BEFORE starting
# - Runs AMRFinderPlus (--plus) on /all_assemblies
# - Recompiles master_AMRFinder_report.tsv robustly
# =========================================================

# -------------------------
# CONFIG
# -------------------------
BASE_DIR="$SEN_ROOT"

# ✅ Use all_assemblies now
ASSEMBLY_DIR="${SEN_AMR_ASSEMBLY_DIR:-$BASE_DIR/all_assemblies}"

OUT_DIR="${SEN_AMRFINDER_OUT:-$BASE_DIR/AMRFinderplus}"
LOG_DIR="$OUT_DIR/logs"

# Threadripper tuning (JOBS * THREADS ≈ 128)
JOBS="${SEN_AMRFINDER_JOBS:-64}"
THREADS="${SEN_AMRFINDER_THREADS:-2}"

ORG="Salmonella"
MIN_FASTA_BYTES=1000

TS="$(date +%F_%H%M%S)"
ARCHIVE_DIR="$OUT_DIR/_archive_${TS}"

MASTER_REPORT="$OUT_DIR/master_AMRFinder_report.tsv"
MANIFEST="$OUT_DIR/assemblies_manifest.tsv"
JOBLOG="$OUT_DIR/amrfinder_joblog_${TS}.tsv"
FAIL_LIST="$OUT_DIR/failed_amrfinder_${TS}.txt"
HITCOUNT_TSV="$OUT_DIR/amrfinder_hitcounts_${TS}.tsv"
MERGE_REPORT="$OUT_DIR/master_build_report_${TS}.tsv"

# -------------------------
# PRECHECKS
# -------------------------
command -v amrfinder >/dev/null 2>&1 || { echo "[ERROR] amrfinder not in PATH"; exit 1; }
command -v parallel >/dev/null 2>&1 || { echo "[ERROR] GNU parallel not in PATH"; exit 1; }
[[ -d "$ASSEMBLY_DIR" ]] || { echo "[ERROR] Missing ASSEMBLY_DIR: $ASSEMBLY_DIR"; exit 1; }

mkdir -p "$OUT_DIR" "$LOG_DIR"
: > "$FAIL_LIST"

echo "[INFO] ASSEMBLY_DIR: $ASSEMBLY_DIR"
echo "[INFO] OUT_DIR:      $OUT_DIR"
echo "[INFO] JOBS:         $JOBS"
echo "[INFO] THREADS/job:  $THREADS"
echo "[INFO] ARCHIVE_DIR:  $ARCHIVE_DIR"

# -------------------------
# ARCHIVE previous outputs (move, don't delete)
# -------------------------
echo "[INFO] Archiving previous outputs..."
mkdir -p "$ARCHIVE_DIR" "$ARCHIVE_DIR/logs"

shopt -s nullglob
for f in "$OUT_DIR"/SRR*.tsv "$OUT_DIR"/*.tsv "$OUT_DIR"/*.txt "$OUT_DIR"/*.log; do
  [[ -e "$f" ]] && mv -f "$f" "$ARCHIVE_DIR"/ || true
done
for f in "$LOG_DIR"/*; do
  [[ -e "$f" ]] && mv -f "$f" "$ARCHIVE_DIR/logs"/ || true
done
shopt -u nullglob

# recreate clean logs dir
rm -rf "$LOG_DIR"
mkdir -p "$LOG_DIR"
: > "$FAIL_LIST"

# -------------------------
# Update DB
# -------------------------
echo "[INFO] Updating AMRFinderPlus database..."
amrfinder -u

# -------------------------
# Build manifest from all_assemblies
# - Sample is SRRxxxx extracted from filename if present
# - Otherwise uses basename without extension
# -------------------------
echo "[INFO] Building assemblies manifest..."
: > "$MANIFEST"

find "$ASSEMBLY_DIR" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) -print0 | \
  awk -v RS='\0' 'NF{print}' | \
  while IFS= read -r fa; do
    [[ -s "$fa" ]] || continue
    [[ "$(stat -c%s "$fa")" -ge "$MIN_FASTA_BYTES" ]] || continue

    bn="$(basename "$fa")"
    # extract SRR id if present
    sample="$(echo "$bn" | grep -Eo 'SRR[0-9]+' | head -n 1 || true)"
    if [[ -z "$sample" ]]; then
      sample="${bn%%.*}"
    fi
    printf "%s\t%s\n" "$sample" "$fa" >> "$MANIFEST"
  done

TOTAL_ASM=$(wc -l < "$MANIFEST" || echo 0)
echo "[INFO] FASTA assemblies discovered: $TOTAL_ASM"
[[ "$TOTAL_ASM" -gt 0 ]] || { echo "[ERROR] No assemblies found in $ASSEMBLY_DIR"; exit 1; }

# -------------------------
# Run AMRFinderPlus
# -------------------------
echo "[INFO] Running AMRFinderPlus on ALL assemblies..."
export OUT_DIR LOG_DIR THREADS ORG FAIL_LIST

# feed fasta paths; derive sample id from manifest first column (safer)
cut -f1,2 "$MANIFEST" | tr '\n' '\0' | \
parallel -0 -j "$JOBS" --colsep '\t' --joblog "$JOBLOG" --eta '
  sample="{1}"
  fa="{2}"
  out="$OUT_DIR/${sample}.tsv"
  log="$LOG_DIR/${sample}.log"

  amrfinder -n "$fa" --organism "$ORG" --plus --threads "$THREADS" -o "$out" >"$log" 2>&1

  if [[ ! -s "$out" ]]; then
    echo "[FAIL] $sample (empty output)" >> "'"$FAIL_LIST"'"
    exit 1
  fi
'

# -------------------------
# Merge master (ROBUST)
# - Only uses *.tsv that are per-sample outputs (from manifest)
# -------------------------
echo "[INFO] Building master report..."

# build list of expected output files from manifest
mapfile -t OUTS < <(cut -f1 "$MANIFEST" | awk '{print "'"$OUT_DIR"'/"$1".tsv"}')

# determine expected header NF by MODE
EXPECTED_NF="$(
  for f in "${OUTS[@]}"; do
    [[ -s "$f" ]] || continue
    awk -F'\t' 'NR==1{print NF; exit}' "$f"
  done | sort | uniq -c | sort -nr | awk 'NR==1{print $2}'
)"
[[ -n "$EXPECTED_NF" ]] || { echo "[ERROR] Could not determine expected header NF."; exit 1; }
echo "[INFO] Merge will use header NF=$EXPECTED_NF"

# pick a header template
HEADER_TSV=""
for f in "${OUTS[@]}"; do
  [[ -s "$f" ]] || continue
  nf="$(awk -F'\t' 'NR==1{print NF; exit}' "$f")"
  if [[ "$nf" == "$EXPECTED_NF" ]]; then
    HEADER_TSV="$f"
    break
  fi
done
[[ -n "$HEADER_TSV" ]] || { echo "[ERROR] No header template TSV found with NF=$EXPECTED_NF"; exit 1; }

TMP_MASTER="$MASTER_REPORT.tmp"
echo -e "Sample\t$(head -n 1 "$HEADER_TSV")" > "$TMP_MASTER"

kept=0
skipped=0

for f in "${OUTS[@]}"; do
  [[ -s "$f" ]] || { ((skipped++)); continue; }
  nf="$(awk -F'\t' 'NR==1{print NF; exit}' "$f")"
  if [[ "$nf" != "$EXPECTED_NF" ]]; then
    ((skipped++))
    continue
  fi
  sample="$(basename "$f" .tsv)"
  tail -n +2 "$f" | awk -v s="$sample" 'BEGIN{OFS="\t"} {print s, $0}' >> "$TMP_MASTER"
  ((kept++))
done

mv -f "$TMP_MASTER" "$MASTER_REPORT"

# hit counts
echo -e "Sample\tN_hits" > "$HITCOUNT_TSV"
for f in "${OUTS[@]}"; do
  [[ -s "$f" ]] || continue
  sample="$(basename "$f" .tsv)"
  n=$(( $(wc -l < "$f") - 1 )); (( n < 0 )) && n=0
  echo -e "${sample}\t${n}" >> "$HITCOUNT_TSV"
done

# merge report
{
  echo -e "metric\tvalue"
  echo -e "expected_header_NF\t$EXPECTED_NF"
  echo -e "master_header_NF_including_Sample\t$((EXPECTED_NF+1))"
  echo -e "assemblies_total\t$TOTAL_ASM"
  echo -e "outputs_kept\t$kept"
  echo -e "outputs_skipped\t$skipped"
  echo -e "master_lines\t$(wc -l < "$MASTER_REPORT")"
} > "$MERGE_REPORT"

echo
echo "✅ DONE"
echo "Archived old run to: $ARCHIVE_DIR"
echo "Manifest:   $MANIFEST"
echo "Master:     $MASTER_REPORT"
echo "Joblog:     $JOBLOG"
echo "Fail list:  $FAIL_LIST"
echo "Hitcounts:  $HITCOUNT_TSV"
echo "Merge rpt:  $MERGE_REPORT"
