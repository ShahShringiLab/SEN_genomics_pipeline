#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
PASS_LIST="$REPO_ROOT/tests/smoke_test/smoke_pass_list.txt"
SNIPPY_ENV="${SEN_SMOKE_SNIPPY_ENV:-sen_snippy}"
GUBBINS_ENV="${SEN_SMOKE_GUBBINS_ENV:-sen_gubbins}"
IQTREE_ENV="${SEN_SMOKE_IQTREE_ENV:-sen_iqtree}"

SNIPPY_YAML="$REPO_ROOT/environments/05_snippy_snpeff.yaml"
GUBBINS_YAML="$REPO_ROOT/environments/06_gubbins.yaml"
IQTREE_YAML="$REPO_ROOT/environments/07_iqtree.yaml"

mkdir -p "$WORKDIR" "$WORKDIR/tmp"
export SEN_ROOT="$WORKDIR"
export TMPDIR="$WORKDIR/tmp"
export SEN_REFERENCE_DIR="$WORKDIR/reference"
export SEN_REFERENCE_FNA="$WORKDIR/reference/reference.fna"
export SEN_REFERENCE_GBK="$WORKDIR/reference/reference.gbk"
export SEN_KRAKEN_CLEAN_DIR="$WORKDIR/Kraken_cleanup/clean_trimmed_fastq"
export SEN_PASS_LIST="$PASS_LIST"
export SEN_SNIPPY_OUT="$WORKDIR/Snippy_output"
export SEN_GUBBINS_OUT="$WORKDIR/gubbin"
export SEN_IQTREE_OUT="$WORKDIR/iqtree_final"

# Beast-safe defaults while remaining overrideable on other machines.
export SEN_SNIPPY_JOBS="${SEN_SNIPPY_JOBS:-4}"
export SEN_SNIPPY_THREADS_PER_JOB="${SEN_SNIPPY_THREADS_PER_JOB:-8}"
export SEN_GUBBINS_THREADS="${SEN_GUBBINS_THREADS:-16}"
export SEN_IQTREE_THREADS="${SEN_IQTREE_THREADS:-16}"
export SEN_IQTREE_MEMORY="${SEN_IQTREE_MEMORY:-32G}"
export SEN_IQTREE_BOOTSTRAPS="${SEN_IQTREE_BOOTSTRAPS:-1000}"

# shellcheck disable=SC1091
source "$REPO_ROOT/tests/smoke_test/checkpoint_lib.sh"

EXPECTED_ISOLATES="$(grep -cve '^[[:space:]]*$' "$PASS_LIST")"
EXPECTED_ALIGNMENT=$((EXPECTED_ISOLATES + 1))

echo "=================================================="
echo " SEN smoke test: reference -> Snippy -> Gubbins -> IQ-TREE"
echo "=================================================="
echo "[INFO] Repo:        $REPO_ROOT"
echo "[INFO] Workdir:     $WORKDIR"
echo "[INFO] Isolates:    $EXPECTED_ISOLATES"
echo "[INFO] + reference: $EXPECTED_ALIGNMENT sequences"

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
    echo "[WARN] $env_name is incomplete; rebuilding."
    conda env remove -n "$env_name" -y
    create_env "$env_name" "$yaml"
  fi
}

validate_reference() {
  [[ -s "$SEN_REFERENCE_FNA" && -s "$SEN_REFERENCE_GBK" ]] || return 1
  grep -q 'NC_011294.1' "$SEN_REFERENCE_FNA" || return 1
  [[ -s "$SEN_REFERENCE_DIR/.sen_reference_verified" ]]
}

validate_snippy() {
  [[ -s "$SEN_SNIPPY_OUT/senbio_core.full.aln" && -s "$SEN_SNIPPY_OUT/senbio_core.vcf" ]] || return 1
  [[ "$(grep -c '^>' "$SEN_SNIPPY_OUT/senbio_core.full.aln" || true)" -eq "$EXPECTED_ALIGNMENT" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$SEN_SNIPPY_OUT/$s/snps.tab" && -s "$SEN_SNIPPY_OUT/$s/snps.vcf" ]] || return 1
  done < "$PASS_LIST"
}

GUBBINS_FILTERED="$SEN_GUBBINS_OUT/senbio_res.filtered_polymorphic_sites.fasta"
validate_gubbins() {
  [[ -s "$GUBBINS_FILTERED" ]] || return 1
  [[ "$(grep -c '^>' "$GUBBINS_FILTERED" || true)" -eq "$EXPECTED_ALIGNMENT" ]]
}

CLEAN_ALN="$SEN_GUBBINS_OUT/clean_final_alignment.fasta"
validate_cleanup() {
  [[ -s "$CLEAN_ALN" ]] || return 1
  [[ "$(grep -c '^>' "$CLEAN_ALN" || true)" -eq "$EXPECTED_ALIGNMENT" ]]
}

TREE="$SEN_IQTREE_OUT/SSLAB_FINAL.treefile"
validate_iqtree() {
  [[ -s "$TREE" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    grep -q "$s" "$TREE" || return 1
  done < "$PASS_LIST"
}

echo
echo "--------------------------------------------------"
echo "[STAGE] 1/5 Reference bootstrap"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "reference_p125109" validate_reference     "$REPO_ROOT/scripts/setup_reference.sh" "$SNIPPY_YAML"; then
  ensure_env "$SNIPPY_ENV" "$SNIPPY_YAML" snippy
  conda run --no-capture-output -n "$SNIPPY_ENV"     env SEN_ROOT="$SEN_ROOT" SEN_REFERENCE_DIR="$SEN_REFERENCE_DIR"         SEN_REFERENCE_FNA="$SEN_REFERENCE_FNA" SEN_REFERENCE_GBK="$SEN_REFERENCE_GBK"     bash "$REPO_ROOT/scripts/setup_reference.sh"
  checkpoint_require_valid "reference_p125109" validate_reference
  checkpoint_mark "reference_p125109" "$REPO_ROOT/scripts/setup_reference.sh" "$SNIPPY_YAML"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 2/5 Snippy + snippy-core"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "snippy_core" validate_snippy     "$REPO_ROOT/workflow/04_core_snp/13_snippy.sh" "$SNIPPY_YAML" "$PASS_LIST"; then
  ensure_env "$SNIPPY_ENV" "$SNIPPY_YAML" snippy
  conda run --no-capture-output -n "$SNIPPY_ENV"     env       SEN_ROOT="$SEN_ROOT" SEN_REFERENCE_DIR="$SEN_REFERENCE_DIR"       SEN_REFERENCE_FNA="$SEN_REFERENCE_FNA" SEN_KRAKEN_CLEAN_DIR="$SEN_KRAKEN_CLEAN_DIR"       SEN_PASS_LIST="$SEN_PASS_LIST" SEN_SNIPPY_OUT="$SEN_SNIPPY_OUT"       SEN_SNIPPY_JOBS="$SEN_SNIPPY_JOBS"       SEN_SNIPPY_THREADS_PER_JOB="$SEN_SNIPPY_THREADS_PER_JOB"       TMPDIR="$TMPDIR"       bash "$REPO_ROOT/workflow/04_core_snp/13_snippy.sh"
  checkpoint_require_valid "snippy_core" validate_snippy
  checkpoint_mark "snippy_core"     "$REPO_ROOT/workflow/04_core_snp/13_snippy.sh" "$SNIPPY_YAML" "$PASS_LIST"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 3/5 Gubbins"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "gubbins" validate_gubbins     "$REPO_ROOT/workflow/05_phylogeny/17_gubbins.sh" "$GUBBINS_YAML"     "$SEN_SNIPPY_OUT/senbio_core.full.aln"; then
  ensure_env "$GUBBINS_ENV" "$GUBBINS_YAML" run_gubbins.py
  conda run --no-capture-output -n "$GUBBINS_ENV"     env       SEN_ROOT="$SEN_ROOT" SEN_CORE_ALIGNMENT="$SEN_SNIPPY_OUT/senbio_core.full.aln"       SEN_GUBBINS_OUT="$SEN_GUBBINS_OUT" SEN_GUBBINS_THREADS="$SEN_GUBBINS_THREADS"       TMPDIR="$TMPDIR"       bash "$REPO_ROOT/workflow/05_phylogeny/17_gubbins.sh"
  checkpoint_require_valid "gubbins" validate_gubbins
  checkpoint_mark "gubbins"     "$REPO_ROOT/workflow/05_phylogeny/17_gubbins.sh" "$GUBBINS_YAML"     "$SEN_SNIPPY_OUT/senbio_core.full.aln"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 4/5 Final-tree cleanup"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "gubbins_cleanup" validate_cleanup     "$REPO_ROOT/workflow/05_phylogeny/18_gubbins_cleanup.sh"     "$REPO_ROOT/config/final_tree_exclusions.txt" "$GUBBINS_FILTERED"; then
  bash "$REPO_ROOT/workflow/05_phylogeny/18_gubbins_cleanup.sh"
  checkpoint_require_valid "gubbins_cleanup" validate_cleanup
  checkpoint_mark "gubbins_cleanup"     "$REPO_ROOT/workflow/05_phylogeny/18_gubbins_cleanup.sh"     "$REPO_ROOT/config/final_tree_exclusions.txt" "$GUBBINS_FILTERED"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 5/5 IQ-TREE"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "iqtree" validate_iqtree     "$REPO_ROOT/workflow/05_phylogeny/19_iqtree.sh" "$IQTREE_YAML" "$CLEAN_ALN"; then
  ensure_env "$IQTREE_ENV" "$IQTREE_YAML" iqtree2
  conda run --no-capture-output -n "$IQTREE_ENV"     env       SEN_ROOT="$SEN_ROOT" SEN_GUBBINS_OUT="$SEN_GUBBINS_OUT"       SEN_IQTREE_INPUT="$CLEAN_ALN" SEN_IQTREE_OUT="$SEN_IQTREE_OUT"       SEN_IQTREE_THREADS="$SEN_IQTREE_THREADS" SEN_IQTREE_MEMORY="$SEN_IQTREE_MEMORY"       SEN_IQTREE_BOOTSTRAPS="$SEN_IQTREE_BOOTSTRAPS" TMPDIR="$TMPDIR"       bash "$REPO_ROOT/workflow/05_phylogeny/19_iqtree.sh"
  checkpoint_require_valid "iqtree" validate_iqtree
  checkpoint_mark "iqtree"     "$REPO_ROOT/workflow/05_phylogeny/19_iqtree.sh" "$IQTREE_YAML" "$CLEAN_ALN"
fi

echo
echo "=================================================="
echo " CORE-SNP / PHYLOGENY SMOKE SUMMARY"
echo "=================================================="
echo "Reference:             NC_011294.1"
echo "Snippy isolate dirs:   $EXPECTED_ISOLATES/$EXPECTED_ISOLATES"
echo "Core alignment taxa:   $(grep -c '^>' "$SEN_SNIPPY_OUT/senbio_core.full.aln")/$EXPECTED_ALIGNMENT"
echo "Gubbins alignment taxa:$(grep -c '^>' "$GUBBINS_FILTERED")/$EXPECTED_ALIGNMENT"
echo "Final alignment taxa:  $(grep -c '^>' "$CLEAN_ALN")/$EXPECTED_ALIGNMENT"
echo "IQ-TREE tree:          $TREE"
echo
echo "SMOKE CORE-SNP/PHYLOGENY RESULT: PASS"
