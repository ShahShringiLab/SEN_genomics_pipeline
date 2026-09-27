#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

WORKDIR="${SEN_SMOKE_ROOT:-$REPO_ROOT/tests/smoke_test/work}"
STATS_ENV="${SEN_SMOKE_STATS_ENV:-sen_stats}"
STATS_YAML="$REPO_ROOT/environments/12_stats.yaml"
PANAROO_ENV="${SEN_SMOKE_PANAROO_ENV:-sen_panaroo}"
PANAROO_YAML="$REPO_ROOT/environments/10_panaroo.yaml"

export SEN_ROOT="$WORKDIR"
export TMPDIR="$WORKDIR/tmp"
export SEN_METADATA_FILE="$REPO_ROOT/metadata/SEN_Genomes.csv"
export SEN_CLADE_METADATA="$WORKDIR/metadata/final_clade_metadata.tsv"
export SEN_COLLECTION_YEAR_ITOL="$WORKDIR/itol_4_collection_year.txt"
export SEN_SOURCE_ITOL="$WORKDIR/itol_1_source.txt"
export SEN_REFERENCE_GBK="$WORKDIR/reference/reference.gbk"
export SEN_FINAL_CONTIGS_DIR="$WORKDIR/Final_Contigs_Only"
export SEN_GUBBINS_OUT="$WORKDIR/gubbin"
export SEN_GUBBINS_RECOMB_GFF="$WORKDIR/gubbin/senbio_res.recombination_predictions.gff"
export SEN_SNP_MASTER_REPORT="$WORKDIR/Snippy_output/SENBIO_RECOVERED_REPORT.csv"
export SEN_SNP_INTEGRATED_ROOT="$WORKDIR/iqtree_final/Snps_clade_integrated"
export SEN_SNP_CLADE_OUT="$SEN_SNP_INTEGRATED_ROOT"
export SEN_SNP_SOURCE_OUT="$WORKDIR/iqtree_final/SNPs_by_source_modest_variation"
export SEN_PANAROO_MATRIX="$WORKDIR/Panaroo_Run/results/gene_presence_absence.csv"
export SEN_PANAROO_GENE_DATA="$WORKDIR/Panaroo_Run/results/gene_data.csv"
export SEN_PANAROO_CLADE_OUT="$WORKDIR/iqtree_final/Panaroo_byclade_robust"

mkdir -p "$WORKDIR/tmp"

# shellcheck disable=SC1091
source "$REPO_ROOT/tests/smoke_test/checkpoint_lib.sh"

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
    echo "[WARN] $env_name incomplete; rebuilding."
    conda env remove -n "$env_name" -y
    create_env "$env_name" "$yaml"
  fi
}

validate_inputs() {
  [[ -s "$WORKDIR/AMRFinderplus/master_AMRFinder_report.tsv" ]] || return 1
  [[ -s "$WORKDIR/abricate_results/resfinder/master_resfinder_report.tsv" ]] || return 1
  [[ -s "$WORKDIR/abricate_results/vfdb/master_vfdb_report.tsv" ]] || return 1
  [[ -s "$WORKDIR/abricate_results/plasmidfinder/master_plasmidfinder_report.tsv" ]] || return 1
  [[ -s "$WORKDIR/Panaroo_Run/results/gene_presence_absence.csv" ]] || return 1
  [[ -s "$WORKDIR/reference/reference.gbk" ]] || return 1
}

validate_clades() {
  [[ -d "$WORKDIR/iqtree_final/clade_summary" ]] || return 1
  [[ -d "$WORKDIR/iqtree_final/AMRFinder_clade" ]] || return 1
  [[ -d "$WORKDIR/iqtree_final/Plasmid_clade" ]] || return 1
  [[ -d "$WORKDIR/iqtree_final/Resfinder_clade" ]] || return 1
  [[ -d "$WORKDIR/iqtree_final/VFDB_clade" ]] || return 1
  find "$WORKDIR/iqtree_final/AMRFinder_clade" -type f -size +0c -print -quit | grep -q .
}

validate_panaroo_clade() {
  [[ -d "$SEN_PANAROO_CLADE_OUT" ]] || return 1
  find "$SEN_PANAROO_CLADE_OUT" -type f -size +0c -print -quit | grep -q .
}

validate_snps() {
  [[ -s "$SEN_SNP_MASTER_REPORT" ]] || return 1
  [[ -d "$SEN_SNP_INTEGRATED_ROOT" ]] || return 1
  find "$SEN_SNP_INTEGRATED_ROOT" -type f -name 'Detailed_Feature_Report.csv' -print -quit | grep -q .
  [[ -d "$SEN_SNP_SOURCE_OUT" ]] || return 1
}

validate_inputs || {
  echo "[ERROR] Upstream smoke outputs incomplete. Run the main smoke pipeline first." >&2
  exit 1
}

ensure_env "$STATS_ENV" "$STATS_YAML" python
ensure_env "$PANAROO_ENV" "$PANAROO_YAML" panaroo

echo "=================================================="
echo " SEN smoke test: downstream analysis layer"
echo "=================================================="

echo
echo "--------------------------------------------------"
echo "[STAGE] 1/5 Build smoke analysis inputs"
echo "--------------------------------------------------"
conda run --no-capture-output -n "$STATS_ENV"   env SEN_SMOKE_ROOT="$WORKDIR"   python "$REPO_ROOT/tests/smoke_test/build_analysis_inputs.py"

# Syntax/load coverage for every source file in analysis/clades and analysis/snps.
echo
echo "--------------------------------------------------"
echo "[STAGE] 2/5 Static load check: analysis/clades + analysis/snps"
echo "--------------------------------------------------"
conda run --no-capture-output -n "$STATS_ENV"   python -m compileall -q "$REPO_ROOT/analysis/clades" "$REPO_ROOT/analysis/snps"
conda run --no-capture-output -n "$STATS_ENV"   Rscript -e 'parse(file="analysis/snps/FGA.R"); parse(file="analysis/snps/FGA2.R"); cat("[PASS] R analysis scripts parse\n")'

echo
echo "--------------------------------------------------"
echo "[STAGE] 3/5 Clade analysis runtime"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "analysis_clades" validate_clades     "$REPO_ROOT/analysis/clades/14_clademetadata.py"     "$REPO_ROOT/analysis/clades/14b_cladesummary.py"     "$REPO_ROOT/analysis/clades/15_amr_byclade.py"     "$REPO_ROOT/analysis/clades/16_plasmid_byclade.py"     "$REPO_ROOT/analysis/clades/17_resfinder_byclade.py"     "$REPO_ROOT/analysis/clades/18_vfdb_byclade.py"     "$STATS_YAML"; then

  conda run --no-capture-output -n "$STATS_ENV"     env SEN_ROOT="$SEN_ROOT"     python "$REPO_ROOT/analysis/clades/14_clademetadata.py"

  for script in     14b_cladesummary.py     15_amr_byclade.py     16_plasmid_byclade.py     17_resfinder_byclade.py     18_vfdb_byclade.py; do
    conda run --no-capture-output -n "$STATS_ENV"       env         SEN_ROOT="$SEN_ROOT"         SEN_METADATA_FILE="$SEN_METADATA_FILE"         SEN_CLADE_METADATA="$SEN_CLADE_METADATA"         SEN_COLLECTION_YEAR_ITOL="$SEN_COLLECTION_YEAR_ITOL"         TMPDIR="$TMPDIR"       python "$REPO_ROOT/analysis/clades/$script"
  done

  checkpoint_require_valid "analysis_clades" validate_clades
  checkpoint_mark "analysis_clades"     "$REPO_ROOT/analysis/clades/14b_cladesummary.py"     "$REPO_ROOT/analysis/clades/15_amr_byclade.py"     "$REPO_ROOT/analysis/clades/16_plasmid_byclade.py"     "$REPO_ROOT/analysis/clades/17_resfinder_byclade.py"     "$REPO_ROOT/analysis/clades/18_vfdb_byclade.py"     "$STATS_YAML"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 4/5 Panaroo clade/accessory analysis runtime"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "analysis_panaroo_clade" validate_panaroo_clade     "$REPO_ROOT/workflow/06_assembly_pangenome/23_panaroo_byclade.py"     "$PANAROO_YAML" "$SEN_CLADE_METADATA"; then
  conda run --no-capture-output -n "$PANAROO_ENV"     env       SEN_ROOT="$SEN_ROOT"       SEN_CLADE_METADATA="$SEN_CLADE_METADATA"       SEN_PANAROO_MATRIX="$SEN_PANAROO_MATRIX"       SEN_PANAROO_GENE_DATA="$SEN_PANAROO_GENE_DATA"       SEN_PANAROO_CLADE_OUT="$SEN_PANAROO_CLADE_OUT"       TMPDIR="$TMPDIR"     python "$REPO_ROOT/workflow/06_assembly_pangenome/23_panaroo_byclade.py"

  checkpoint_require_valid "analysis_panaroo_clade" validate_panaroo_clade
  checkpoint_mark "analysis_panaroo_clade"     "$REPO_ROOT/workflow/06_assembly_pangenome/23_panaroo_byclade.py"     "$PANAROO_YAML" "$SEN_CLADE_METADATA"
fi

echo
echo "--------------------------------------------------"
echo "[STAGE] 5/5 SNP analysis runtime"
echo "--------------------------------------------------"
if ! checkpoint_should_skip "analysis_snps" validate_snps     "$REPO_ROOT/analysis/snps/19b_snps_byclade.py"     "$REPO_ROOT/analysis/snps/19c_snps_proteinseq_fetch.py"     "$REPO_ROOT/analysis/snps/19d_snps_unmatched_query_rescue.py"     "$REPO_ROOT/analysis/snps/19e_snps_fga.py"     "$REPO_ROOT/analysis/snps/19f_snps_bysource.py"     "$STATS_YAML"; then

  for script in     19b_snps_byclade.py     19c_snps_proteinseq_fetch.py     19d_snps_unmatched_query_rescue.py     19e_snps_fga.py     19f_snps_bysource.py; do
    conda run --no-capture-output -n "$STATS_ENV"       env         SEN_ROOT="$SEN_ROOT"         SEN_CLADE_METADATA="$SEN_CLADE_METADATA"         SEN_REFERENCE_GBK="$SEN_REFERENCE_GBK"         SEN_FINAL_CONTIGS_DIR="$SEN_FINAL_CONTIGS_DIR"         SEN_GUBBINS_OUT="$SEN_GUBBINS_OUT"         SEN_GUBBINS_RECOMB_GFF="$SEN_GUBBINS_RECOMB_GFF"         SEN_SNP_MASTER_REPORT="$SEN_SNP_MASTER_REPORT"         SEN_SNP_INTEGRATED_ROOT="$SEN_SNP_INTEGRATED_ROOT"         SEN_SNP_CLADE_OUT="$SEN_SNP_CLADE_OUT"         SEN_SOURCE_ITOL="$SEN_SOURCE_ITOL"         SEN_SNP_SOURCE_OUT="$SEN_SNP_SOURCE_OUT"         TMPDIR="$TMPDIR"       python "$REPO_ROOT/analysis/snps/$script"
  done

  checkpoint_require_valid "analysis_snps" validate_snps
  checkpoint_mark "analysis_snps"     "$REPO_ROOT/analysis/snps/19b_snps_byclade.py"     "$REPO_ROOT/analysis/snps/19c_snps_proteinseq_fetch.py"     "$REPO_ROOT/analysis/snps/19d_snps_unmatched_query_rescue.py"     "$REPO_ROOT/analysis/snps/19e_snps_fga.py"     "$REPO_ROOT/analysis/snps/19f_snps_bysource.py"     "$STATS_YAML"
fi

echo
echo "=================================================="
echo " DOWNSTREAM ANALYSIS SMOKE SUMMARY"
echo "=================================================="
echo "analysis/clades: runtime PASS"
echo "Panaroo-by-clade: runtime PASS"
echo "analysis/snps Python chain: runtime PASS"
echo "analysis/snps R scripts: parse PASS"
echo
echo "[NOTE] Four-isolate outputs are runtime/contract tests only; statistical"
echo "[NOTE] significance and biological enrichment are not interpreted."
echo
echo "SMOKE DOWNSTREAM ANALYSIS RESULT: PASS"
