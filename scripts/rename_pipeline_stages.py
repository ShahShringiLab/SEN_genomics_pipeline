#!/usr/bin/env python3
"""Stage-prefixed naming migration for SEN pipeline scripts.

Usage:
  python scripts/rename_pipeline_stages.py --check
  python scripts/rename_pipeline_stages.py --apply

--apply renames tracked scripts with git mv and updates tracked text references.
Run from a clean Git checkout, inspect git diff and execute validators before commit.
"""
from __future__ import annotations
import argparse
import re
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RENAMES={
 "analysis/part_a/00_parse_master_33831.R":"analysis/part_a/A2A_parse_master_33831.R",
 "analysis/part_a/01_phenotype_trends.R":"analysis/part_a/A1_phenotype_trends.R",
 "analysis/part_a/02_build_final_cohort.R":"analysis/part_a/A2B_build_final_cohort.R",
 "analysis/part_a/03_genotype_marker_summary.R":"analysis/part_a/A2C_genotype_marker_summary.R",
 "analysis/part_a/04_source_profiles.R":"analysis/part_a/A2D_source_profiles.R",
 "analysis/part_a/05_diversity_composition.R":"analysis/part_a/A2E_diversity_composition.R",
 "analysis/part_a/06_master_33831_analysis.R":"analysis/part_a/A2F_master_33831_analysis.R",
 "analysis/part_a/07_master_full_diversity.R":"analysis/part_a/A2G_master_full_diversity.R",
 "analysis/part_a/07_selection_qc_weights_table.R":"analysis/part_a/A2H_selection_qc_weights_table.R",
 "analysis/part_a/08_plot_full_diversity_original_style.R":"analysis/part_a/A2I_plot_full_diversity_original_style.R",
 "analysis/part_a/08_table_ss1_post2018_genotype.R":"analysis/part_a/A2J_table_ss1_post2018_genotype.R",
 "analysis/part_a/09_temporal_top5_and_enriched_metadata.R":"analysis/part_a/A2K_temporal_top5_and_enriched_metadata.R",
 "analysis/part_a/10_international_animal_source_plot.R":"analysis/part_a/A2L_international_animal_source_plot.R",
 "workflow/01_reads/01_download_fastq.sh":"workflow/01_reads/B1A_download_fastq.sh",
 "workflow/02_qc/02_fastqc_raw.sh":"workflow/02_qc/B1B_fastqc_raw.sh",
 "workflow/02_qc/03_fastqc_trimmed.sh":"workflow/02_qc/B1C_fastqc_trimmed.sh",
 "workflow/02_qc/04_estimate_coverage.sh":"workflow/02_qc/B1D_estimate_coverage.sh",
 "workflow/02_qc/05_two_step_filter.sh":"workflow/02_qc/B1E_two_step_filter.sh",
 "workflow/02_qc/06_kraken_cleanup.sh":"workflow/02_qc/B1F_kraken_cleanup.sh",
 "workflow/02_qc/07_presnippy_qc.py":"workflow/02_qc/B1G_presnippy_qc.py",
 "workflow/02_qc/08_qc_filter.py":"workflow/02_qc/B1H_qc_filter.py",
 "workflow/03_typing/09_trim_fastp.sh":"workflow/03_typing/B1I_trim_fastp.sh",
 "workflow/03_typing/10_seqsero2.sh":"workflow/03_typing/B2A_seqsero2.sh",
 "workflow/03_typing/11_sistr.sh":"workflow/03_typing/B2B_sistr.sh",
 "workflow/03_typing/12_mlst.sh":"workflow/03_typing/B2C_mlst.sh",
 "workflow/04_core_snp/13_snippy.sh":"workflow/04_core_snp/B3A_snippy.sh",
 "workflow/04_core_snp/14_snpeff.sh":"workflow/04_core_snp/B3B_snpeff.sh",
 "workflow/04_core_snp/15_snpeff_parse.py":"workflow/04_core_snp/B3C_snpeff_parse.py",
 "workflow/04_core_snp/16_genome_audit.py":"workflow/04_core_snp/B3D_genome_audit.py",
 "workflow/05_phylogeny/17_gubbins.sh":"workflow/05_phylogeny/B3E_gubbins.sh",
 "workflow/05_phylogeny/18_gubbins_cleanup.sh":"workflow/05_phylogeny/B3F_gubbins_cleanup.sh",
 "workflow/05_phylogeny/19_iqtree.sh":"workflow/05_phylogeny/B3G_iqtree.sh",
 "workflow/05_phylogeny/20_tree_clades.sh":"workflow/05_phylogeny/B3H_tree_clades.sh",
 "workflow/06_assembly_pangenome/21_shovill.sh":"workflow/06_assembly_pangenome/B4A_shovill.sh",
 "workflow/06_assembly_pangenome/22_prokka.py":"workflow/06_assembly_pangenome/B4B_prokka.py",
 "workflow/06_assembly_pangenome/23_panaroo.sh":"workflow/06_assembly_pangenome/B4C_panaroo.sh",
 "workflow/06_assembly_pangenome/23_panaroo_byclade.py":"workflow/06_assembly_pangenome/B4D_panaroo_byclade.py",
 "workflow/07_amr_vf_plasmid/24_amrfinderplus.sh":"workflow/07_amr_vf_plasmid/B5A_amrfinderplus.sh",
 "workflow/07_amr_vf_plasmid/25_abricate.sh":"workflow/07_amr_vf_plasmid/B5B_abricate.sh",
 "analysis/clades/14_clademetadata.py":"analysis/clades/B3I_clademetadata.py",
 "analysis/clades/14b_cladesummary.py":"analysis/clades/B3J_cladesummary.py",
 "analysis/clades/15_amr_byclade.py":"analysis/clades/B5C_amr_byclade.py",
 "analysis/clades/16_plasmid_byclade.py":"analysis/clades/B5D_plasmid_byclade.py",
 "analysis/clades/17_resfinder_byclade.py":"analysis/clades/B5E_resfinder_byclade.py",
 "analysis/clades/18_vfdb_byclade.py":"analysis/clades/B5F_vfdb_byclade.py",
 "analysis/snps/19b_snps_byclade.py":"analysis/snps/B3K_snps_byclade.py",
 "analysis/snps/19c_snps_proteinseq_fetch.py":"analysis/snps/B3L_snps_proteinseq_fetch.py",
 "analysis/snps/19d_snps_unmatched_query_rescue.py":"analysis/snps/B3M_snps_unmatched_query_rescue.py",
 "analysis/snps/19e_snps_fga.py":"analysis/snps/B3N_snps_fga.py",
 "analysis/snps/19f_snps_bysource.py":"analysis/snps/B3O_snps_bysource.py",
}
def git(*args):
    return subprocess.check_output(["git",*args],cwd=ROOT,text=True).strip()

def main():
    ap=argparse.ArgumentParser()
    grp=ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--check",action="store_true")
    grp.add_argument("--apply",action="store_true")
    args=ap.parse_args()
    tracked=git("ls-files","-z").split("\0")
    tracked=[x for x in tracked if x]
    missing=[x for x in RENAMES if x not in tracked]
    collisions=[v for v in RENAMES.values() if v in tracked]
    if missing or collisions: raise SystemExit(f"Invalid mapping: missing={missing}; collisions={collisions}")
    print(f"[CHECK] {len(RENAMES)} file renames; tracked references in {len(tracked)} files")
    if args.check: return
    if git("status","--porcelain"):
        raise SystemExit("[STOP] Git worktree must be clean; commit/stash changes first")
    # Move all scripts. No content or output modifications to archived files.
    for old,new in RENAMES.items():
        subprocess.check_call(["git","mv",old,new],cwd=ROOT)
    # Update references across the tracked text files (avoid data and preserved provenance).
    ref_files=[]
    for old in tracked:
        current=RENAMES.get(old,old)
        p=ROOT/current
        if not p.is_file() or p.suffix.lower() in {".csv",".tsv",".xlsx",".png",".pdf",".fa",".fasta",".gz",".rmd"}:
            continue
        if current.startswith(("data/","metadata/","provenance/")): continue
        raw=p.read_bytes()
        if b"\0" in raw: continue
        text=raw.decode("utf-8")
        rewritten=text
        for before,after in RENAMES.items():
            rewritten=rewritten.replace(before,after)
        # Bare basenames commonly occur in runners and test scripts.
        for before,after in RENAMES.items():
            rewritten=rewritten.replace(Path(before).name,Path(after).name)
        if rewritten!=text:
            p.write_text(rewritten,encoding="utf-8")
            ref_files.append(current)
    print(f"[PASS] Renamed {len(RENAMES)} scripts and updated {len(ref_files)} text files.")
    print("[NEXT] git diff --check; run tests/validation; inspect diff; commit explicitly.")
if __name__=="__main__":main()
