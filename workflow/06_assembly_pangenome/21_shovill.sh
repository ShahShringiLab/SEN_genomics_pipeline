#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd shovill
require_cmd parallel

INPUT_DIR="${SEN_KRAKEN_CLEAN_DIR:-$SEN_ROOT/Kraken_cleanup/clean_trimmed_fastq}"
OUTPUT_DIR="${SEN_SHOVILL_OUT:-$SEN_ROOT/shovill_assemblies}"
FINAL_CONTIGS="${SEN_FINAL_CONTIGS_DIR:-$SEN_ROOT/Final_Contigs_Only}"
MANIFEST="${SEN_SHOVILL_MANIFEST:-$SEN_ROOT/shovill_manifest.tsv}"
LOG_DIR="${SEN_SHOVILL_LOG_DIR:-$SEN_ROOT/logs/shovill}"
TMP_ROOT="${SEN_SHOVILL_TMP:-$SEN_ROOT/tmp_shovill}"

JOBS="${SEN_SHOVILL_JOBS:-4}"
CPUS_PER_JOB="${SEN_SHOVILL_CPUS_PER_JOB:-8}"
RAM_PER_JOB_GB="${SEN_SHOVILL_RAM_PER_JOB_GB:-16}"
MINLEN="${SEN_SHOVILL_MINLEN:-200}"
DEPTH="${SEN_SHOVILL_DEPTH:-0}"
ASSEMBLER="${SEN_SHOVILL_ASSEMBLER:-spades}"

require_dir "$INPUT_DIR"
mkdir -p "$OUTPUT_DIR" "$FINAL_CONTIGS" "$LOG_DIR" "$TMP_ROOT"

echo "[INFO] Shovill input:       $INPUT_DIR"
echo "[INFO] Shovill output:      $OUTPUT_DIR"
echo "[INFO] Final contigs:       $FINAL_CONTIGS"
echo "[INFO] jobs=$JOBS cpus/job=$CPUS_PER_JOB ram/job=${RAM_PER_JOB_GB}G"
echo "[INFO] assembler=$ASSEMBLER minlen=$MINLEN depth=$DEPTH"

: > "$MANIFEST"
shopt -s nullglob
for r1 in "$INPUT_DIR"/*_pure_1.fastq.gz; do
  sample="$(basename "$r1" _pure_1.fastq.gz)"
  r2="$INPUT_DIR/${sample}_pure_2.fastq.gz"
  [[ -s "$r2" ]] || {
    echo "[ERROR] Missing R2 for $sample: $r2" >&2
    exit 1
  }
  printf '%s\t%s\t%s\n' "$sample" "$r1" "$r2" >> "$MANIFEST"
done

[[ -s "$MANIFEST" ]] || {
  echo "[ERROR] No paired Kraken-cleaned FASTQs found in $INPUT_DIR" >&2
  exit 1
}

assemble_one() {
  local sample="$1" r1="$2" r2="$3"
  local outdir="$OUTPUT_DIR/$sample"
  local tmpdir="$TMP_ROOT/$sample"
  local final="$FINAL_CONTIGS/${sample}.fasta"
  local outlog="$LOG_DIR/${sample}.out.log"
  local errlog="$LOG_DIR/${sample}.err.log"

  if [[ -s "$outdir/contigs.fa" ]]; then
    ln -sf "$(realpath "$outdir/contigs.fa")" "$final"
    echo "[SKIP] Shovill $sample"
    return 0
  fi

  rm -rf "$outdir" "$tmpdir"
  mkdir -p "$tmpdir"

  echo "[RUN] Shovill $sample"
  if ! shovill       --outdir "$outdir"       --R1 "$r1" --R2 "$r2"       --cpus "$CPUS_PER_JOB"       --ram "$RAM_PER_JOB_GB"       --assembler "$ASSEMBLER"       --minlen "$MINLEN"       --depth "$DEPTH"       --tmpdir "$tmpdir"       --force       >"$outlog" 2>"$errlog"; then
    rm -rf "$outdir"
    echo "[ERROR] Shovill failed for $sample; see $errlog" >&2
    return 1
  fi

  [[ -s "$outdir/contigs.fa" ]] || {
    echo "[ERROR] Shovill produced no contigs for $sample" >&2
    return 1
  }

  ln -sf "$(realpath "$outdir/contigs.fa")" "$final"
  rm -rf "$tmpdir"
  echo "[DONE] Shovill $sample"
}

export -f assemble_one
export OUTPUT_DIR FINAL_CONTIGS LOG_DIR TMP_ROOT CPUS_PER_JOB RAM_PER_JOB_GB
export MINLEN DEPTH ASSEMBLER

parallel --colsep '\t' -j "$JOBS" --halt now,fail=1   assemble_one {1} {2} {3} :::: "$MANIFEST"

expected="$(wc -l < "$MANIFEST")"
observed="$(find "$FINAL_CONTIGS" -maxdepth 1 -type l -o -type f | grep -E '\.(fasta|fa|fna)$' | wc -l || true)"

echo "[INFO] Shovill assemblies expected: $expected"
echo "[INFO] Final contig files present:  $observed"
echo "[INFO] Shovill stage complete."
