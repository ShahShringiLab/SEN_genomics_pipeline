#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
SRR_LIST="$REPO_ROOT/tests/smoke_test/smoke_srrs.txt"

SHOVILL_ENV="${SEN_SMOKE_SHOVILL_ENV:-sen_shovill}"
PROKKA_ENV="${SEN_SMOKE_PROKKA_ENV:-sen_prokka}"
PANAROO_ENV="${SEN_SMOKE_PANAROO_ENV:-sen_panaroo}"

SHOVILL_YAML="$REPO_ROOT/environments/08_shovill.yaml"
PROKKA_YAML="$REPO_ROOT/environments/09_prokka.yaml"
PANAROO_YAML="$REPO_ROOT/environments/10_panaroo.yaml"

mkdir -p "$WORKDIR" "$WORKDIR/tmp"
export SEN_ROOT="$WORKDIR"
export TMPDIR="$WORKDIR/tmp"
export SEN_KRAKEN_CLEAN_DIR="$WORKDIR/Kraken_cleanup/clean_trimmed_fastq"
export SEN_SHOVILL_OUT="$WORKDIR/shovill_assemblies"
export SEN_FINAL_CONTIGS_DIR="$WORKDIR/Final_Contigs_Only"
export SEN_PROKKA_OUT="$WORKDIR/prokka_annotations"
export SEN_PANAROO_OUT="$WORKDIR/Panaroo_Run/results"

export SEN_SHOVILL_JOBS="${SEN_SHOVILL_JOBS:-4}"
export SEN_SHOVILL_CPUS_PER_JOB="${SEN_SHOVILL_CPUS_PER_JOB:-8}"
export SEN_SHOVILL_RAM_PER_JOB_GB="${SEN_SHOVILL_RAM_PER_JOB_GB:-16}"
export SEN_PROKKA_JOBS="${SEN_PROKKA_JOBS:-4}"
export SEN_PROKKA_CPUS_PER_JOB="${SEN_PROKKA_CPUS_PER_JOB:-4}"
export SEN_PANAROO_THREADS="${SEN_PANAROO_THREADS:-16}"

# shellcheck disable=SC1091
source "$REPO_ROOT/tests/smoke_test/checkpoint_lib.sh"

EXPECTED="$(grep -cve '^[[:space:]]*$' "$SRR_LIST")"

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda is not available in PATH." >&2
  exit 1
}

create_env() {
  local env_name="$1" yaml="$2"
  echo "[INFO] Creating $env_name from $yaml ..."
  if conda env create --help 2>&1 | grep -q -- '--solver'; then
    conda env create --solver=libmamba -f "$yaml"
  else
    conda env create -f "$yaml"
  fi
}

ensure_env() {
  local env_name="$1" yaml="$2" cmd="$3"
  if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
    create_env "$env_name" "$yaml"
  elif ! conda run -n "$env_name" command -v "$cmd" >/dev/null 2>&1; then
    echo "[WARN] $env_name exists but $cmd is unavailable; rebuilding."
    conda env remove -n "$env_name" -y
    create_env "$env_name" "$yaml"
  fi
}

validate_shovill() {
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$SEN_SHOVILL_OUT/$s/contigs.fa" ]] || return 1
    [[ -s "$SEN_FINAL_CONTIGS_DIR/$s.fasta" ]] || return 1
  done < "$SRR_LIST"
}

validate_prokka() {
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$SEN_PROKKA_OUT/$s.gff" ]] || return 1
    grep -q '^##gff-version' "$SEN_PROKKA_OUT/$s.gff" || return 1
  done < "$SRR_LIST"
}

PAN_GPA="$SEN_PANAROO_OUT/gene_presence_absence.csv"
PAN_CORE="$SEN_PANAROO_OUT/core_gene_alignment.aln"
validate_panaroo() {
  [[ -s "$PAN_GPA" && -s "$PAN_CORE" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    grep -q "$s" "$PAN_GPA" || return 1
  done < "$SRR_LIST"
}

echo "=================================================="
echo " SEN smoke test: Shovill -> Prokka -> Panaroo"
echo "=================================================="
echo "[INFO] Workdir: $WORKDIR"
echo "[INFO] Samples: $EXPECTED"

echo
echo "--------------------------------------------------"
echo "[STAGE] 1/3 Shovill assemblies"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "shovill" validate_shovill     "$REPO_ROOT/workflow/06_assembly_pangenome/21_shovill.sh" "$SHOVILL_YAML" "$SRR_LIST"; then
  ensure_env "$SHOVILL_ENV" "$SHOVILL_YAML" shovill
  conda run --no-capture-output -n "$SHOVILL_ENV"     env       SEN_ROOT="$SEN_ROOT" SEN_KRAKEN_CLEAN_DIR="$SEN_KRAKEN_CLEAN_DIR"       SEN_SHOVILL_OUT="$SEN_SHOVILL_OUT" SEN_FINAL_CONTIGS_DIR="$SEN_FINAL_CONTIGS_DIR"       SEN_SHOVILL_JOBS="$SEN_SHOVILL_JOBS"       SEN_SHOVILL_CPUS_PER_JOB="$SEN_SHOVILL_CPUS_PER_JOB"       SEN_SHOVILL_RAM_PER_JOB_GB="$SEN_SHOVILL_RAM_PER_JOB_GB"       TMPDIR="$TMPDIR"       bash "$REPO_ROOT/workflow/06_assembly_pangenome/21_shovill.sh"
  checkpoint_require_valid "shovill" validate_shovill
  checkpoint_mark "shovill"     "$REPO_ROOT/workflow/06_assembly_pangenome/21_shovill.sh" "$SHOVILL_YAML" "$SRR_LIST"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 2/3 Prokka annotation"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "prokka" validate_prokka     "$REPO_ROOT/workflow/06_assembly_pangenome/22_prokka.py" "$PROKKA_YAML" "$SRR_LIST"; then
  ensure_env "$PROKKA_ENV" "$PROKKA_YAML" prokka
  conda run --no-capture-output -n "$PROKKA_ENV"     env       SEN_ROOT="$SEN_ROOT" SEN_FINAL_CONTIGS_DIR="$SEN_FINAL_CONTIGS_DIR"       SEN_PROKKA_OUT="$SEN_PROKKA_OUT"       SEN_PROKKA_JOBS="$SEN_PROKKA_JOBS"       SEN_PROKKA_CPUS_PER_JOB="$SEN_PROKKA_CPUS_PER_JOB"       TMPDIR="$TMPDIR"       python "$REPO_ROOT/workflow/06_assembly_pangenome/22_prokka.py"
  checkpoint_require_valid "prokka" validate_prokka
  checkpoint_mark "prokka"     "$REPO_ROOT/workflow/06_assembly_pangenome/22_prokka.py" "$PROKKA_YAML" "$SRR_LIST"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 3/3 Panaroo"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "panaroo" validate_panaroo     "$REPO_ROOT/workflow/06_assembly_pangenome/23_panaroo.sh" "$PANAROO_YAML" "$SRR_LIST"; then
  ensure_env "$PANAROO_ENV" "$PANAROO_YAML" panaroo

  # Panaroo 1.6.0 parses Prokka-derived temporary FASTA with SeqIO format
  # "fasta". Biopython 1.87 made leading-comment handling an error, so reject
  # environments that resolve that incompatible parser behavior.
  if ! conda run -n "$PANAROO_ENV" python -c '
import Bio, sys
assert sys.version_info[:2] == (3, 10)
assert Bio.__version__ == "1.86"
' >/dev/null 2>&1; then
    echo "[WARN] $PANAROO_ENV has an incompatible Python/Biopython runtime; rebuilding."
    conda env remove -n "$PANAROO_ENV" -y
    create_env "$PANAROO_ENV" "$PANAROO_YAML"
  fi

  rm -rf "$SEN_PANAROO_OUT"
  conda run --no-capture-output -n "$PANAROO_ENV"     env       SEN_ROOT="$SEN_ROOT" SEN_PROKKA_OUT="$SEN_PROKKA_OUT"       SEN_PANAROO_OUT="$SEN_PANAROO_OUT" SEN_PANAROO_THREADS="$SEN_PANAROO_THREADS"       TMPDIR="$TMPDIR"       bash "$REPO_ROOT/workflow/06_assembly_pangenome/23_panaroo.sh"
  checkpoint_require_valid "panaroo" validate_panaroo
  checkpoint_mark "panaroo"     "$REPO_ROOT/workflow/06_assembly_pangenome/23_panaroo.sh" "$PANAROO_YAML" "$SRR_LIST"
fi

echo
echo "=================================================="
echo " ASSEMBLY / PANGENOME SMOKE SUMMARY"
echo "=================================================="
echo "Shovill assemblies: $EXPECTED/$EXPECTED"
echo "Prokka GFFs:        $EXPECTED/$EXPECTED"
echo "Panaroo matrix:     $PAN_GPA"
echo "Panaroo core aln:   $PAN_CORE"
echo
echo "SMOKE ASSEMBLY/PANGENOME RESULT: PASS"
