#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
SRR_LIST="$REPO_ROOT/tests/smoke_test/smoke_srrs.txt"
ENV_NAME="${SEN_SMOKE_AMR_ENV:-sen_amr_abricate}"
ENV_YAML="$REPO_ROOT/environments/11_amr_abricate.yaml"

mkdir -p "$WORKDIR" "$WORKDIR/tmp"
export SEN_ROOT="$WORKDIR"
export TMPDIR="$WORKDIR/tmp"
export SEN_FINAL_CONTIGS_DIR="$WORKDIR/Final_Contigs_Only"
export SEN_AMR_DB_ROOT="$WORKDIR/databases"
export SEN_AMRFINDER_DB_ROOT="$WORKDIR/databases/amrfinderplus"
export SEN_AMR_PROVENANCE_DIR="$WORKDIR/databases/provenance"
export SEN_AMRFINDER_OUT="$WORKDIR/AMRFinderplus"
export SEN_ABRICATE_OUT="$WORKDIR/abricate_results"

export SEN_AMRFINDER_JOBS="${SEN_AMRFINDER_JOBS:-4}"
export SEN_AMRFINDER_THREADS="${SEN_AMRFINDER_THREADS:-4}"
export SEN_ABRICATE_JOBS="${SEN_ABRICATE_JOBS:-4}"
export SEN_ABRICATE_THREADS="${SEN_ABRICATE_THREADS:-2}"
export SEN_ABRICATE_MINCOV="${SEN_ABRICATE_MINCOV:-80}"
export SEN_ABRICATE_MINID="${SEN_ABRICATE_MINID:-90}"

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
  local env_name="$1" yaml="$2"
  if ! conda env list | awk '{print $1}' | grep -qx "$env_name"; then
    create_env "$env_name" "$yaml"
  elif ! conda run -n "$env_name" command -v amrfinder >/dev/null 2>&1 ||        ! conda run -n "$env_name" command -v abricate >/dev/null 2>&1; then
    echo "[WARN] $env_name is incomplete; rebuilding."
    conda env remove -n "$env_name" -y
    create_env "$env_name" "$yaml"
  fi
}

validate_inputs() {
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$SEN_FINAL_CONTIGS_DIR/$s.fasta" ]] || return 1
  done < "$SRR_LIST"
}

validate_databases() {
  [[ -e "$SEN_AMRFINDER_DB_ROOT/latest" ]] || return 1
  local db
  db="$(readlink -f "$SEN_AMRFINDER_DB_ROOT/latest" 2>/dev/null || true)"
  [[ -n "$db" && -d "$db" && -s "$db/version.txt" ]] || return 1
  [[ -s "$SEN_AMR_PROVENANCE_DIR/amrfinderplus.tsv" ]] || return 1
  [[ -s "$SEN_AMR_PROVENANCE_DIR/abricate_databases.tsv" ]] || return 1
  for x in resfinder vfdb plasmidfinder; do
    awk -v db="$x" 'BEGIN{FS="\t"} NR>1 && $1==db && $2+0>0 {ok=1} END{exit !ok}'       "$SEN_AMR_PROVENANCE_DIR/abricate_databases.tsv" || return 1
  done
}

validate_amrfinder() {
  [[ -s "$SEN_AMRFINDER_OUT/master_AMRFinder_report.tsv" ]] || return 1
  [[ -s "$SEN_AMRFINDER_OUT/database_version.txt" ]] || return 1
  local s
  while IFS= read -r s; do
    [[ -n "$s" ]] || continue
    [[ -s "$SEN_AMRFINDER_OUT/$s.tsv" ]] || return 1
  done < "$SRR_LIST"
}

validate_abricate() {
  [[ -s "$SEN_ABRICATE_OUT/database_inventory.tsv" ]] || return 1
  local db s
  for db in resfinder vfdb plasmidfinder; do
    [[ -s "$SEN_ABRICATE_OUT/$db/master_${db}_report.tsv" ]] || return 1
    while IFS= read -r s; do
      [[ -n "$s" ]] || continue
      [[ -s "$SEN_ABRICATE_OUT/$db/${s}_${db}.tsv" ]] || return 1
    done < "$SRR_LIST"
  done
}

validate_inputs || {
  echo "[ERROR] Final Shovill contigs are incomplete. Run assembly/pangenome smoke first." >&2
  exit 1
}

echo "=================================================="
echo " SEN smoke test: AMRFinderPlus + ABRicate"
echo "=================================================="
echo "[INFO] Workdir: $WORKDIR"
echo "[INFO] Samples: $EXPECTED"
echo "[INFO] ABRicate thresholds: coverage >=$SEN_ABRICATE_MINCOV%, identity >=$SEN_ABRICATE_MINID%"

ensure_env "$ENV_NAME" "$ENV_YAML"

echo
echo "--------------------------------------------------"
echo "[STAGE] 1/3 AMR database bootstrap/provenance"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "amr_databases" validate_databases     "$REPO_ROOT/scripts/setup_amr_databases.sh" "$ENV_YAML"; then
  conda run --no-capture-output -n "$ENV_NAME"     env       SEN_ROOT="$SEN_ROOT" SEN_AMR_DB_ROOT="$SEN_AMR_DB_ROOT"       SEN_AMRFINDER_DB_ROOT="$SEN_AMRFINDER_DB_ROOT"       SEN_AMR_PROVENANCE_DIR="$SEN_AMR_PROVENANCE_DIR"       TMPDIR="$TMPDIR"       bash "$REPO_ROOT/scripts/setup_amr_databases.sh"
  checkpoint_require_valid "amr_databases" validate_databases
  checkpoint_mark "amr_databases" "$REPO_ROOT/scripts/setup_amr_databases.sh" "$ENV_YAML"
fi

export SEN_AMRFINDER_DB="$(readlink -f "$SEN_AMRFINDER_DB_ROOT/latest")"
echo "[INFO] Frozen AMRFinderPlus DB: $SEN_AMRFINDER_DB"

echo
echo "--------------------------------------------------"
echo "[STAGE] 2/3 AMRFinderPlus"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "amrfinderplus" validate_amrfinder     "$REPO_ROOT/workflow/07_amr_vf_plasmid/24_amrfinderplus.sh" "$ENV_YAML"     "$SEN_AMR_PROVENANCE_DIR/amrfinderplus.tsv" "$SRR_LIST"; then
  conda run --no-capture-output -n "$ENV_NAME"     env       SEN_ROOT="$SEN_ROOT" SEN_FINAL_CONTIGS_DIR="$SEN_FINAL_CONTIGS_DIR"       SEN_AMR_ASSEMBLY_DIR="$SEN_FINAL_CONTIGS_DIR"       SEN_AMRFINDER_OUT="$SEN_AMRFINDER_OUT" SEN_AMRFINDER_DB="$SEN_AMRFINDER_DB"       SEN_AMRFINDER_JOBS="$SEN_AMRFINDER_JOBS"       SEN_AMRFINDER_THREADS="$SEN_AMRFINDER_THREADS"       TMPDIR="$TMPDIR"       bash "$REPO_ROOT/workflow/07_amr_vf_plasmid/24_amrfinderplus.sh"
  checkpoint_require_valid "amrfinderplus" validate_amrfinder
  checkpoint_mark "amrfinderplus"     "$REPO_ROOT/workflow/07_amr_vf_plasmid/24_amrfinderplus.sh" "$ENV_YAML"     "$SEN_AMR_PROVENANCE_DIR/amrfinderplus.tsv" "$SRR_LIST"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 3/3 ABRicate ResFinder/VFDB/PlasmidFinder"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "abricate" validate_abricate     "$REPO_ROOT/workflow/07_amr_vf_plasmid/25_abricate.sh" "$ENV_YAML"     "$SEN_AMR_PROVENANCE_DIR/abricate_databases.tsv" "$SRR_LIST"; then
  conda run --no-capture-output -n "$ENV_NAME"     env       SEN_ROOT="$SEN_ROOT" SEN_FINAL_CONTIGS_DIR="$SEN_FINAL_CONTIGS_DIR"       SEN_ABRICATE_ASSEMBLY_DIR="$SEN_FINAL_CONTIGS_DIR"       SEN_ABRICATE_OUT="$SEN_ABRICATE_OUT"       SEN_ABRICATE_JOBS="$SEN_ABRICATE_JOBS"       SEN_ABRICATE_THREADS="$SEN_ABRICATE_THREADS"       SEN_ABRICATE_MINCOV="$SEN_ABRICATE_MINCOV"       SEN_ABRICATE_MINID="$SEN_ABRICATE_MINID"       TMPDIR="$TMPDIR"       bash "$REPO_ROOT/workflow/07_amr_vf_plasmid/25_abricate.sh"
  checkpoint_require_valid "abricate" validate_abricate
  checkpoint_mark "abricate"     "$REPO_ROOT/workflow/07_amr_vf_plasmid/25_abricate.sh" "$ENV_YAML"     "$SEN_AMR_PROVENANCE_DIR/abricate_databases.tsv" "$SRR_LIST"
fi

echo
echo "=================================================="
echo " AMR / VIRULENCE / PLASMID SMOKE SUMMARY"
echo "=================================================="
echo "AMRFinderPlus reports: $EXPECTED/$EXPECTED"
echo "AMRFinder database:"
cat "$SEN_AMRFINDER_OUT/database_version.txt" | sed 's/^/  /'
echo
echo "ABRicate database inventory:"
awk 'NR==1 || $1=="resfinder" || $1=="vfdb" || $1=="plasmidfinder"'   "$SEN_ABRICATE_OUT/database_inventory.tsv"
echo
echo "ABRicate thresholds: coverage >=$SEN_ABRICATE_MINCOV%, identity >=$SEN_ABRICATE_MINID%"
echo
echo "SMOKE AMR/VF/PLASMID RESULT: PASS"
