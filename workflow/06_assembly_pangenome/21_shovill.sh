#!/bin/bash
set -euo pipefail

# ==============================================================================
# ULTRAFAST + STABLE Massive Salmonella Assembly (PURE READS; Kraken-cleaned)
# Uses your previous directory layout:
#   INPUT_DIR=/home/samuelajulo/SENBio/Final/Kraken_cleanup/clean_trimmed_fastq
#   OUTPUT_DIR=shovill_assemblies
#   FINAL_CONTIGS=Final_Contigs_Only
# Optimized scheduling: 24 jobs x 5 threads (120 threads) + safety checks
# Temp: prefers /dev/shm if enough free space, else falls back to NVMe/local
# ==============================================================================

# ---------------------------
# 0) CONFIG (MATCHES YOUR PREVIOUS WORKING PATHS)
# ---------------------------
INPUT_DIR="/home/samuelajulo/SENBio/Final/Kraken_cleanup/clean_trimmed_fastq"
OUTPUT_DIR="shovill_assemblies"
FINAL_CONTIGS="Final_Contigs_Only"
MANIFEST="samples_pure.tsv"

# Throughput tuning
CONCURRENT_JOBS=24
CORES_PER_JOB=5
RAM_PER_JOB_GB=20

# Assembly knobs
MINLEN=200
DEPTH=0
ASSEMBLER="spades"          # fallback: skesa
SPADES_OPTS=""

# Paths
WD="$(pwd)"
RUN_ID="$(date +%F_%H%M%S)"
LOG_DIR="$WD/logs"
JOBLOG="$WD/parallel_joblog_${RUN_ID}.tsv"
FAIL_LIST="$WD/failed_samples_${RUN_ID}.txt"

# Temp strategy: use /dev/shm if sufficiently free, else NVMe/local
TMP_SHM="/dev/shm/tmp_shovill_${RUN_ID}"
TMP_NVME="$WD/tmp_shovill_${RUN_ID}"
MIN_SHM_FREE_GB=200
TMP_ROOT="$TMP_NVME"  # default, may switch to shm if safe

# System tweaks
ulimit -n 20480
export OMP_NUM_THREADS=1

# ---------------------------
# 1) PRE-FLIGHT CHECKS
# ---------------------------
echo "--- PRE-FLIGHT CHECKS ---"
need_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "[ERROR] Missing in PATH: $1"; exit 1; }; }

need_cmd parallel
need_cmd shovill
need_cmd java
if [[ "$ASSEMBLER" == "spades" ]]; then
  need_cmd spades.py
fi

mkdir -p "$OUTPUT_DIR" "$FINAL_CONTIGS" "$LOG_DIR" "$TMP_NVME"

echo "[INFO] Run ID:   $RUN_ID"
echo "[INFO] INPUT:    $INPUT_DIR"
echo "[INFO] OUTPUT:   $OUTPUT_DIR"
echo "[INFO] CONTIGS:  $FINAL_CONTIGS"
echo "[INFO] JOBLOG:   $JOBLOG"
echo "[INFO] Jobs: $CONCURRENT_JOBS | Threads/job: $CORES_PER_JOB | RAM/job: ${RAM_PER_JOB_GB}G"
echo

echo "--- SHOVILL TOOLCHAIN CHECK ---"
shovill --check
echo

# Input sanity
[[ -d "$INPUT_DIR" ]] || { echo "[ERROR] INPUT_DIR not found: $INPUT_DIR"; exit 1; }
N_R1=$(ls -1 "$INPUT_DIR"/*_pure_1.fastq.gz 2>/dev/null | wc -l || true)
[[ "$N_R1" -gt 0 ]] || { echo "[ERROR] No *_pure_1.fastq.gz files found in $INPUT_DIR"; exit 1; }

# Choose TMP_ROOT based on /dev/shm free space
if df -BG /dev/shm >/dev/null 2>&1; then
  SHM_FREE_GB=$(df -BG /dev/shm | awk 'NR==2{gsub(/G/,"",$4); print $4}')
  if [[ "${SHM_FREE_GB:-0}" -ge "$MIN_SHM_FREE_GB" ]]; then
    TMP_ROOT="$TMP_SHM"
    mkdir -p "$TMP_ROOT"
    echo "[INFO] Using /dev/shm tmp: $TMP_ROOT (free ${SHM_FREE_GB}G)"
  else
    TMP_ROOT="$TMP_NVME"
    echo "[WARN] /dev/shm free ${SHM_FREE_GB}G < ${MIN_SHM_FREE_GB}G → using NVMe tmp: $TMP_ROOT"
  fi
else
  TMP_ROOT="$TMP_NVME"
  echo "[WARN] Could not read /dev/shm → using NVMe tmp: $TMP_ROOT"
fi
echo

# ---------------------------
# 2) BUILD MANIFEST (PURE READS)
# ---------------------------
echo "--- BUILDING MANIFEST FROM KRAKEN-CLEANED READS ---"
: > "$MANIFEST"
shopt -s nullglob

paired=0
missing=0

for r1 in "$INPUT_DIR"/*_pure_1.fastq.gz; do
  sample=$(basename "$r1" _pure_1.fastq.gz)
  r2="$INPUT_DIR/${sample}_pure_2.fastq.gz"

  if [[ -s "$r2" ]]; then
    printf "%s\t%s\t%s\n" "$sample" "$r1" "$r2" >> "$MANIFEST"
    paired=$((paired+1))
  else
    echo "[WARN] Missing R2 for $sample (expected: $r2) — skipping"
    missing=$((missing+1))
  fi
done

[[ "$paired" -gt 0 ]] || { echo "[ERROR] No paired pure FASTQs found in $INPUT_DIR"; exit 1; }
echo "[INFO] Samples to process: $paired | Missing R2 skipped: $missing"
echo

# ---------------------------
# 3) SMOKE TEST (FIRST SAMPLE)
# ---------------------------
echo "--- SMOKE TEST (first sample) ---"
first_sample=$(awk 'NR==1{print $1}' "$MANIFEST")
first_r1=$(awk 'NR==1{print $2}' "$MANIFEST")
first_r2=$(awk 'NR==1{print $3}' "$MANIFEST")

shovill --outdir "$OUTPUT_DIR/_SMOKE_${first_sample}" \
  --R1 "$first_r1" --R2 "$first_r2" \
  --cpus "$CORES_PER_JOB" --ram "$RAM_PER_JOB_GB" \
  --assembler "$ASSEMBLER" --minlen "$MINLEN" --depth "$DEPTH" \
  --tmpdir "$TMP_ROOT/_SMOKE_${first_sample}" --force \
  > "$LOG_DIR/_SMOKE_${first_sample}.out.log" 2> "$LOG_DIR/_SMOKE_${first_sample}.err.log" || {
    echo "[ERROR] Smoke test failed. See: $LOG_DIR/_SMOKE_${first_sample}.err.log" >&2
    exit 1
  }

if [[ ! -s "$OUTPUT_DIR/_SMOKE_${first_sample}/contigs.fa" ]]; then
  echo "[ERROR] Smoke test produced no contigs. See: $LOG_DIR/_SMOKE_${first_sample}.err.log" >&2
  exit 1
fi
echo "[SMOKE] OK: contigs created for $first_sample"
echo

# ---------------------------
# 4) ASSEMBLY FUNCTION
# ---------------------------
assemble_one() {
  sample="$1"
  r1="$2"
  r2="$3"

  outdir="$OUTPUT_DIR/$sample"
  tmpdir="$TMP_ROOT/$sample"
  outlog="$LOG_DIR/${sample}.shovill.out.log"
  errlog="$LOG_DIR/${sample}.shovill.err.log"

  mkdir -p "$outdir" "$tmpdir"

  # Resume
  if [[ -s "$outdir/contigs.fa" ]]; then
    ln -sf "$(realpath "$outdir/contigs.fa")" "$FINAL_CONTIGS/${sample}.fasta" 2>/dev/null || true
    echo "[SKIP] $sample"
    return 0
  fi

  cmd=(shovill
    --outdir "$outdir"
    --R1 "$r1" --R2 "$r2"
    --cpus "$CORES_PER_JOB"
    --ram "$RAM_PER_JOB_GB"
    --assembler "$ASSEMBLER"
    --minlen "$MINLEN"
    --depth "$DEPTH"
    --tmpdir "$tmpdir"
    --force
  )

  [[ -n "$SPADES_OPTS" ]] && cmd+=(--opts "$SPADES_OPTS")

  if "${cmd[@]}" > "$outlog" 2> "$errlog"; then
    if [[ -s "$outdir/contigs.fa" ]]; then
      ln -sf "$(realpath "$outdir/contigs.fa")" "$FINAL_CONTIGS/${sample}.fasta" 2>/dev/null || true
      rm -rf "$tmpdir" 2>/dev/null || true
      echo "[OK] $sample"
      return 0
    fi
  fi

  rm -rf "$tmpdir" 2>/dev/null || true
  echo "[FAIL] $sample"
  return 1
}

export -f assemble_one
export OUTPUT_DIR FINAL_CONTIGS TMP_ROOT LOG_DIR CORES_PER_JOB RAM_PER_JOB_GB MINLEN ASSEMBLER SPADES_OPTS DEPTH

# ---------------------------
# 5) EXECUTION
# ---------------------------
echo "--- RUNNING ASSEMBLIES ($CONCURRENT_JOBS jobs x $CORES_PER_JOB threads) ---"
parallel --colsep '\t' \
  -j "$CONCURRENT_JOBS" \
  --joblog "$JOBLOG" \
  --resume-failed \
  --eta \
  assemble_one {1} {2} {3} \
  :::: "$MANIFEST" \
  || true

echo "--- PIPELINE FINISHED ---"

# ---------------------------
# 6) SUMMARY
# ---------------------------
: > "$FAIL_LIST"
fails=0
ok=0

while IFS=$'\t' read -r sample r1 r2; do
  if [[ -s "$OUTPUT_DIR/$sample/contigs.fa" ]]; then
    ok=$((ok+1))
  else
    echo "$sample" >> "$FAIL_LIST"
    fails=$((fails+1))
  fi
done < "$MANIFEST"

echo "[INFO] Successful contigs: $ok"
echo "[INFO] Failures: $fails (see $FAIL_LIST)"
echo "[INFO] Final contigs symlinks: $FINAL_CONTIGS/"
echo "[INFO] Logs: $LOG_DIR/"
echo "[INFO] Joblog: $JOBLOG"
echo "[INFO] TMP_ROOT used: $TMP_ROOT"

