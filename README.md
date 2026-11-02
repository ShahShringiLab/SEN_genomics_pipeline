# Salmonella Enteritidis Genomics Pipeline

Reproducible publication workflow for the *Salmonella enterica* serovar Enteritidis (SEN) analyses described in the associated manuscript.

The repository implements the complete study workflow in two connected parts:

- **Part A** — phenotype resistance trends, QR/FQR genotype epidemiology, source/origin classification, marker diversity, and construction/validation of the final 3,434-genome cohort.
- **Part B** — read retrieval and QC, Enteritidis confirmation, core-SNP phylogenomics, recombination filtering, clade analysis, assembly/pangenome analysis, and AMR/virulence/plasmid analyses.

The original exploratory Part A notebook is preserved under provenance/. Linux-executable Part A and Part B scripts are organized together under workflow/.

## Repository structure

    SEN_genomics_pipeline/
    ├── README.md
    ├── CITATION.cff
    ├── analysis/
    │   ├── part_a/
    │   ├── clades/
    │   └── snps/
    ├── config/
    ├── data/
    │   └── part_a/
    ├── docs/
    ├── environments/
    ├── metadata/
    ├── provenance/
    ├── reference/
    ├── scripts/
    ├── tests/
    └── workflow/

## Study workflow

### Part A1 — phenotype resistance trends

**Existing script filenames (unchanged):**

- `workflow/A1_Phenotype/01_phenotype_trends.R`


Input:

    data/part_a/phenotype/FDASEN_AMRPhenotype2.csv

The production analysis reproduces the Human-versus-Chicken temporal analysis of ciprofloxacin resistance, nalidixic-acid reduced susceptibility/resistance, and the combined QR/FQR phenotype metric with 95% confidence intervals.

Run:

    conda run --no-capture-output -n sen_stats \
      env SEN_ROOT="$PWD" \
      Rscript workflow/A1_Phenotype/01_phenotype_trends.R

### Part A2 — genotype epidemiology and cohort construction

**Existing script filenames (unchanged):**

- `workflow/A2_Genotype/00_parse_master_33831.R`
- `workflow/A2_Genotype/02_build_final_cohort.R`
- `workflow/A2_Genotype/03_genotype_marker_summary.R`
- `workflow/A2_Genotype/04_source_profiles.R`
- `workflow/A2_Genotype/05_diversity_composition.R`
- `workflow/A2_Genotype/06_master_33831_analysis.R`
- `workflow/A2_Genotype/07_master_full_diversity.R`
- `workflow/A2_Genotype/07_selection_qc_weights_table.R`
- `workflow/A2_Genotype/08_plot_full_diversity_original_style.R`
- `workflow/A2_Genotype/08_table_ss1_post2018_genotype.R`
- `workflow/A2_Genotype/09_temporal_top5_and_enriched_metadata.R`
- `workflow/A2_Genotype/10_international_animal_source_plot.R`


Inputs:

    data/part_a/genotype/US-SEN_12192025.csv
    data/part_a/selection/Selected_updated_Isolate.csv
    data/part_a/selection/Finalupdated_SEN_Genomes1.csv

Part A2:

1. derives collection year and broad QR/FQR genotype status from the master genotype/metadata table;
2. reproduces the Human/Chicken/Other and US/non-US source logic;
3. applies the final 3,434-isolate selection key;
4. validates rebuilt QR/FQR status, collection year and final source against the committed final cohort;
5. summarizes QR/FQR marker frequency, marker burden, temporal/source patterns, and multi-marker combinations;
6. reproduces final source/origin profiles, including International and Other-US animal-source composition;
7. quantifies QR/FQR diversity and composition using Hill numbers, bootstrap confidence intervals, and permutation tests.

Run the complete Part A workflow with:

    bash scripts/run_part_a.sh

A successful run ends with:

    SEN PART A: PASS

That PASS is emitted only after both the canonical-input/handoff validator and
the generated-output validator succeed. The latter verifies the rebuilt
3,434-genome cohort, the expected three QR-status audit differences, final
source counts, phenotype output dimensions, diversity source coverage, and key
figure artifacts.

Generated outputs are written beneath results/part_a/ and are intentionally not version-controlled.

The archived master genotype CSV is preserved unchanged. Its known CSV quoting/escaped-comma artifacts are normalized only in memory by the Part A reader before strict row/column validation.

The broad QR/FQR marker rule reconstructed from the original notebook differs from the committed final QR_genotype for three isolates. The workflow writes these discrepancies to an audit table and uses the committed final QR_genotype for downstream production analyses rather than silently changing either record.

Fresh-Linux instructions: docs/PART_A_FRESH_LINUX.md.

### Part B1 — read retrieval, trimming and quality control

**Existing script filenames (unchanged):**

- `workflow/B1_Reads_QC/01_download_fastq.sh`
- `workflow/B1_Reads_QC/02_fastqc_raw.sh`
- `workflow/B1_Reads_QC/03_fastqc_trimmed.sh`
- `workflow/B1_Reads_QC/04_estimate_coverage.sh`
- `workflow/B1_Reads_QC/05_two_step_filter.sh`
- `workflow/B1_Reads_QC/06_kraken_cleanup.sh`
- `workflow/B1_Reads_QC/07_presnippy_qc.py`
- `workflow/B1_Reads_QC/08_qc_filter.py`
- `workflow/B1_Reads_QC/09_trim_fastp.sh`


The final Part A cohort provides the 3,434 selected accessions entering Part B.

1. Download Illumina reads with SRA Toolkit.
2. Trim and quality-filter paired-end reads with fastp.
3. Summarize raw and trimmed read quality with FastQC and MultiQC.
4. Classify reads with Kraken2 and retain Enterobacteriaceae reads (TaxID 543).
5. Require post-filter sequencing depth >=30x.

### Part B2 — species confirmation and typing

**Existing script filenames (unchanged):**

- `workflow/B2_Typing/10_seqsero2.sh`
- `workflow/B2_Typing/11_sistr.sh`
- `workflow/B2_Typing/12_mlst.sh`


1. SeqSero2 k-mer serotyping.
2. Rapid SKESA assembly followed by SISTR.
3. Retain genomes confirmed as *Salmonella* Enteritidis by both workflows.
4. MLST typing.

The final confirmed study cohort contains 3,306 study isolates.

### Part B3 — core-genome SNP and phylogenetic analysis

**Existing script filenames (unchanged):**

- `workflow/B3_Phylogenomics/13_snippy.sh`
- `workflow/B3_Phylogenomics/14_snpeff.sh`
- `workflow/B3_Phylogenomics/15_snpeff_parse.py`
- `workflow/B3_Phylogenomics/16_genome_audit.py`
- `workflow/B3_Phylogenomics/17_gubbins.sh`
- `workflow/B3_Phylogenomics/18_gubbins_cleanup.sh`
- `workflow/B3_Phylogenomics/19_iqtree.sh`
- `workflow/B3_Phylogenomics/20_tree_clades.sh`
- `workflow/B3_Phylogenomics/14_clademetadata.py`
- `workflow/B3_Phylogenomics/14b_cladesummary.py`
- `workflow/B3_Phylogenomics/19b_snps_byclade.py`
- `workflow/B3_Phylogenomics/19c_snps_proteinseq_fetch.py`
- `workflow/B3_Phylogenomics/19d_snps_unmatched_query_rescue.py`
- `workflow/B3_Phylogenomics/19e_snps_fga.py`
- `workflow/B3_Phylogenomics/19f_snps_bysource.py`


1. Map reads to SEN reference P125109 (NC_011294.1) with Snippy.
2. Generate the core-genome alignment with snippy-core.
3. Annotate variants with SnpEff.
4. Mask recombination with Gubbins.
5. Infer the maximum-likelihood phylogeny with IQ-TREE.
6. Validate the final 3,306 study isolates plus P125109 against canonical clade metadata.
7. Perform clade-enriched and clade-defining SNP analyses.

### Part B4 — assembly, annotation and pangenome analysis

**Existing script filenames (unchanged):**

- `workflow/B4_Pangenome/21_shovill.sh`
- `workflow/B4_Pangenome/22_prokka.py`
- `workflow/B4_Pangenome/23_panaroo.sh`
- `workflow/B4_Pangenome/23_panaroo_byclade.py`


1. Generate de novo draft assemblies with Shovill/SPAdes.
2. Annotate assemblies with Prokka.
3. Reconstruct the pangenome with Panaroo in strict mode using MAFFT.
4. Test accessory-gene clade enrichment.
5. Define clade-associated accessory genes as >=90% within clade and <=5% outside clade.

### Part B5 — AMR, virulence and plasmid analysis

**Existing script filenames (unchanged):**

- `workflow/B5_AMR_Virulence_Plasmids/24_amrfinderplus.sh`
- `workflow/B5_AMR_Virulence_Plasmids/25_abricate.sh`
- `workflow/B5_AMR_Virulence_Plasmids/15_amr_byclade.py`
- `workflow/B5_AMR_Virulence_Plasmids/16_plasmid_byclade.py`
- `workflow/B5_AMR_Virulence_Plasmids/17_resfinder_byclade.py`
- `workflow/B5_AMR_Virulence_Plasmids/18_vfdb_byclade.py`


1. AMRFinderPlus with --plus --organism Salmonella.
2. ResFinder, VFDB and PlasmidFinder through ABRicate.
3. Binary presence requires >=80% coverage and >=90% identity.
4. Clade-stratified statistical analyses use omnibus testing, sparse-count handling, pairwise Fisher tests and multiple-testing correction.

## Fresh-Linux execution

Create all declared software environments:

    bash scripts/setup_all_envs.sh

Run Part A:

    bash scripts/run_part_a.sh

Run the Part B smoke test on a new machine before manuscript-scale execution:

    bash tests/smoke_test/run_smoke_pipeline.sh

The core end-to-end smoke workflow has been exercised through QC/typing, core-SNP/phylogeny, assembly/pangenome, and AMR/ABRicate stages. Downstream analysis scripts remain subject to their explicit publication-contract checks.

For the complete manuscript-scale Part B run, see:

    docs/FULL_RUN_FRESH_LINUX.md

## Reproducible environments

Software environments are defined under environments/, with analysis-critical versions pinned for reproducible execution.

Part A uses sen_stats, including R, tidyverse, lubridate and patchwork. Part B uses stage-specific environments for reads/QC, Kraken2, serotyping, MLST, Snippy/SnpEff, Gubbins, IQ-TREE, Shovill, Prokka, Panaroo and AMR screening.

## Canonical data and metadata

Part A final cohort:

    data/part_a/selection/Finalupdated_SEN_Genomes1.csv

Part B canonical starting metadata:

    metadata/SEN_Genomes.csv

Canonical final clade metadata:

    metadata/final_clade_metadata.tsv

The final phylogenomic composition is:

    Clade 1A       169
    Clade 1B       759
    Clade 2        894
    Clade 3      1,357
    Unassigned     127
    Reference        1
    -----------------
    Total         3,307 tips

## Configuration

Scientific thresholds are encoded explicitly in the workflow and analysis scripts and audited by scripts/validate_analysis_contracts.py.

Frozen external database releases are recorded in config/database_sources.env.

Machine-specific paths belong in an untracked config/paths.env, using config/paths.example.env as the template.

## Validation

Run the static publication audit with:

    bash scripts/publication_audit.sh

Part A performs two runtime validation layers: it first validates all committed
inputs and the 3,434-accession Part A -> Part B handoff, then validates generated
Part A outputs before reporting PASS. Part B includes a representative smoke
workflow and canonical metadata/clade-count checks.

## Provenance

The original working genomics repository is ajulojays/SENBio-Final.

The original Part A exploratory analysis notebook is retained as:

    provenance/SEN_Genomes2.Rmd

The notebook is preserved for analytical provenance. The production workflow should be run from the modular scripts under analysis/part_a/, not by executing the exploratory notebook top-to-bottom.

## Release status

This repository contains the complete manuscript-facing SEN workflow implementation spanning Part A and Part B.

A fresh run uses the workflow code, canonical data/metadata, and pinned software/database definitions stored in this release.
