#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
KRAKEN_ENV="${SEN_SMOKE_KRAKEN_ENV:-sen_kraken}"
SEQSERO_ENV="${SEN_SMOKE_SEQSERO_ENV:-sen_seqsero2}"
SISTR_ENV="${SEN_SMOKE_SISTR_ENV:-sen_sistr}"
MLST_ENV="${SEN_SMOKE_MLST_ENV:-sen_mlst}"

KRAKEN_YAML="$REPO_ROOT/environments/02_kraken.yaml"
SEQSERO_YAML="$REPO_ROOT/environments/03_seqsero2.yaml"
SISTR_YAML="$REPO_ROOT/environments/03_sistr.yaml"
MLST_YAML="$REPO_ROOT/environments/04_mlst.yaml"

mkdir -p "$WORKDIR" "$WORKDIR/tmp" "$WORKDIR/tmp/kraken"

export SEN_ROOT="$WORKDIR"
export TMPDIR="$WORKDIR/tmp"
export SEN_KRAKEN_TMP="$WORKDIR/tmp/kraken"
export SEN_KRAKEN_CLEAN_DIR="$WORKDIR/Kraken_cleanup/clean_trimmed_fastq"

# Laptop-safe defaults; override before launch if desired.
export SEN_KRAKEN_PARALLEL_JOBS="${SEN_KRAKEN_PARALLEL_JOBS:-2}"
export SEN_KRAKEN_THREADS="${SEN_KRAKEN_THREADS:-4}"
export SEN_PIGZ_THREADS="${SEN_PIGZ_THREADS:-2}"
export SEN_SEQSERO2_JOBS="${SEN_SEQSERO2_JOBS:-2}"
export SEN_SISTR_JOBS="${SEN_SISTR_JOBS:-2}"
export SEN_SKESA_CORES="${SEN_SKESA_CORES:-4}"
export SEN_SKESA_MEMORY_GB="${SEN_SKESA_MEMORY_GB:-4}"
export SEN_MLST_JOBS="${SEN_MLST_JOBS:-4}"

echo "=================================================="
echo " SEN smoke test: Kraken -> SeqSero2/SISTR -> MLST"
echo "=================================================="
echo "[INFO] Repo:       $REPO_ROOT"
echo "[INFO] Workdir:    $WORKDIR"
echo "[INFO] TMPDIR:     $TMPDIR"

command -v conda >/dev/null 2>&1 || {
  echo "[ERROR] conda is not available in PATH." >&2
  exit 1
}

create_env() {
  local env_name="$1"
  local yaml="$2"

  echo "[INFO] Creating $env_name from $yaml ..."
  if conda env create --help 2>&1 | grep -q -- '--solver'; then
    conda env create --solver=libmamba -f "$yaml"
  else
    echo "[ERROR] This Conda installation does not expose --solver=libmamba." >&2
    exit 1
  fi
}

ensure_env() {
  local env_name="$1"
  local yaml="$2"
  local required_cmd="$3"

  if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
    create_env "$env_name" "$yaml"
  fi

  if ! conda run -n "$env_name" command -v "$required_cmd" >/dev/null 2>&1; then
    echo "[WARN] $env_name exists but $required_cmd is unavailable."
    echo "[INFO] Rebuilding dedicated smoke environment $env_name ..."
    conda env remove -n "$env_name" -y
    create_env "$env_name" "$yaml"
  fi
}

ensure_env "$KRAKEN_ENV" "$KRAKEN_YAML" kraken2
ensure_env "$SEQSERO_ENV" "$SEQSERO_YAML" SeqSero2_package.py
ensure_env "$SISTR_ENV" "$SISTR_YAML" skesa
if ! conda run -n "$SISTR_ENV" command -v sistr >/dev/null 2>&1; then
  echo "[WARN] $SISTR_ENV is missing sistr; rebuilding environment."
  conda env remove -n "$SISTR_ENV" -y
  create_env "$SISTR_ENV" "$SISTR_YAML"
fi
ensure_env "$MLST_ENV" "$MLST_YAML" mlst

TRIM_DIR="$WORKDIR/trimmed_fastq"
if [[ ! -d "$TRIM_DIR" ]] || [[ "$(find "$TRIM_DIR" -maxdepth 1 -name '*_trimmed_1.fastq.gz' | wc -l)" -ne 4 ]]; then
  echo "[ERROR] Expected four trimmed smoke-test read pairs in $TRIM_DIR." >&2
  echo "[ERROR] Run: bash tests/smoke_test/run_reads_qc.sh" >&2
  exit 1
fi

# Reuse an explicitly configured database when valid; otherwise bootstrap the
# pinned publication-facing snapshot automatically.
if [[ -z "${SEN_KRAKEN_DB:-}" && -f "$REPO_ROOT/config/paths.env" ]]; then
  # shellcheck disable=SC1091
  source "$REPO_ROOT/config/paths.env"
fi

DEFAULT_KRAKEN_DB="$WORKDIR/databases/kraken2/k2_standard_20260626"

if [[ -z "${SEN_KRAKEN_DB:-}" ||       ! -s "${SEN_KRAKEN_DB:-}/hash.k2d" ||       ! -s "${SEN_KRAKEN_DB:-}/opts.k2d" ||       ! -s "${SEN_KRAKEN_DB:-}/taxo.k2d" ]]; then
  echo "[INFO] No complete Kraken2 database configured."
  echo "[INFO] Bootstrapping the pinned database snapshot automatically..."

  conda run --no-capture-output -n "$KRAKEN_ENV"     env SEN_ROOT="$SEN_ROOT" TMPDIR="$TMPDIR"     bash "$REPO_ROOT/scripts/setup_kraken_db.sh"

  export SEN_KRAKEN_DB="$DEFAULT_KRAKEN_DB"
fi

if [[ ! -s "$SEN_KRAKEN_DB/hash.k2d" ||       ! -s "$SEN_KRAKEN_DB/opts.k2d" ||       ! -s "$SEN_KRAKEN_DB/taxo.k2d" ]]; then
  echo "[ERROR] Kraken2 database bootstrap did not produce a complete database." >&2
  exit 1
fi

echo "[INFO] Kraken DB:   $SEN_KRAKEN_DB"

echo
echo "--------------------------------------------------"
echo "[STAGE] 1/4 Kraken2 filtering + Enterobacteriaceae extraction"
echo "--------------------------------------------------"
conda run --no-capture-output -n "$KRAKEN_ENV"   env     SEN_ROOT="$SEN_ROOT"     SEN_KRAKEN_DB="$SEN_KRAKEN_DB"     SEN_KRAKEN_TMP="$SEN_KRAKEN_TMP"     SEN_KRAKEN_PARALLEL_JOBS="$SEN_KRAKEN_PARALLEL_JOBS"     SEN_KRAKEN_THREADS="$SEN_KRAKEN_THREADS"     SEN_PIGZ_THREADS="$SEN_PIGZ_THREADS"     TMPDIR="$TMPDIR"     bash "$REPO_ROOT/workflow/02_qc/06_kraken_cleanup.sh"

KRAKEN_LOG="$WORKDIR/Kraken_cleanup/full_db_results.csv"
[[ -s "$KRAKEN_LOG" ]] || { echo "[FAIL] Missing Kraken summary: $KRAKEN_LOG" >&2; exit 1; }

post30="$(awk -F',' 'NR>1 && $3>=30 {n++} END{print n+0}' "$KRAKEN_LOG")"
if [[ "$post30" -ne 4 ]]; then
  echo "[FAIL] Only $post30/4 smoke isolates retain >=30x post-Kraken depth." >&2
  cat "$KRAKEN_LOG"
  exit 1
fi

pure_pairs="$(find "$SEN_KRAKEN_CLEAN_DIR" -maxdepth 1 -name '*_pure_1.fastq.gz' | wc -l)"
[[ "$pure_pairs" -eq 4 ]] || { echo "[FAIL] Expected 4 Kraken-cleaned read pairs, found $pure_pairs." >&2; exit 1; }

echo
echo "--------------------------------------------------"
echo "[STAGE] 2/4 SeqSero2"
echo "--------------------------------------------------"
conda run --no-capture-output -n "$SEROTYPE_ENV"   env     SEN_ROOT="$SEN_ROOT"     SEN_SEQSERO2_JOBS="$SEN_SEQSERO2_JOBS"     TMPDIR="$TMPDIR"     bash "$REPO_ROOT/workflow/03_typing/10_seqsero2.sh"

SEQSERO_SUMMARY="$WORKDIR/seqsero2_results/SeqSero2_summary.tsv"
[[ -s "$SEQSERO_SUMMARY" ]] || { echo "[FAIL] Missing SeqSero2 summary." >&2; exit 1; }

echo
echo "--------------------------------------------------"
echo "[STAGE] 3/4 SKESA + SISTR"
echo "--------------------------------------------------"
conda run --no-capture-output -n "$SEROTYPE_ENV"   env     SEN_ROOT="$SEN_ROOT"     SEN_KRAKEN_CLEAN_DIR="$SEN_KRAKEN_CLEAN_DIR"     SEN_SISTR_JOBS="$SEN_SISTR_JOBS"     SEN_SKESA_CORES="$SEN_SKESA_CORES"     SEN_SKESA_MEMORY_GB="$SEN_SKESA_MEMORY_GB"     TMPDIR="$TMPDIR"     bash "$REPO_ROOT/workflow/03_typing/11_sistr.sh"

SISTR_SUMMARY="$WORKDIR/sistr_results_run/sistr_master_summary.csv"
[[ -s "$SISTR_SUMMARY" ]] || { echo "[FAIL] Missing SISTR summary." >&2; exit 1; }

assembly_count="$(find "$WORKDIR/sistr_results_run/mini_assemblies" -maxdepth 1 -name '*.fasta' | wc -l)"
[[ "$assembly_count" -eq 4 ]] || { echo "[FAIL] Expected 4 SKESA assemblies, found $assembly_count." >&2; exit 1; }

echo
echo "--------------------------------------------------"
echo "[STAGE] 4/4 MLST"
echo "--------------------------------------------------"
conda run --no-capture-output -n "$MLST_ENV"   env     SEN_ROOT="$SEN_ROOT"     SEN_MLST_JOBS="$SEN_MLST_JOBS"     TMPDIR="$TMPDIR"     bash "$REPO_ROOT/workflow/03_typing/12_mlst.sh"

MLST_SUMMARY="$WORKDIR/mlst_results_run/mlst_master_report.csv"
[[ -s "$MLST_SUMMARY" ]] || { echo "[FAIL] Missing MLST summary." >&2; exit 1; }

seqsero_enteritidis="$(grep -ic 'Enteritidis' "$SEQSERO_SUMMARY" || true)"
sistr_enteritidis="$(grep -ic 'Enteritidis' "$SISTR_SUMMARY" || true)"
mlst_rows="$(( $(wc -l < "$MLST_SUMMARY") - 1 ))"

echo
echo "=================================================="
echo " KRAKEN/TYPING SMOKE SUMMARY"
echo "=================================================="
cat "$KRAKEN_LOG"
echo
echo "Kraken-cleaned read pairs: $pure_pairs/4"
echo "Post-Kraken depth >=30x:   $post30/4"
echo "SeqSero2 Enteritidis hits: $seqsero_enteritidis"
echo "SISTR Enteritidis hits:    $sistr_enteritidis"
echo "SKESA assemblies:          $assembly_count/4"
echo "MLST result rows:           $mlst_rows"

if [[ "$seqsero_enteritidis" -lt 4 ]]; then
  echo "[FAIL] SeqSero2 did not report Enteritidis for all four smoke isolates." >&2
  exit 1
fi
if [[ "$sistr_enteritidis" -lt 4 ]]; then
  echo "[FAIL] SISTR did not report Enteritidis for all four smoke isolates." >&2
  exit 1
fi
if [[ "$mlst_rows" -lt 4 ]]; then
  echo "[FAIL] MLST returned fewer than four results." >&2
  exit 1
fi

echo
echo "SMOKE KRAKEN/TYPING RESULT: PASS"
