# Salmonella Enteritidis Genomics Pipeline

Reproducible publication workflow for the *Salmonella enterica* serovar Enteritidis (SEN) analyses described in the associated manuscript.

The repository implements the complete study workflow in two connected parts:

- **Part A** — phenotype resistance trends, QR/FQR genotype epidemiology, source/origin classification, marker diversity, and construction/validation of the final 3,434-genome cohort.
- **Part B** — read retrieval and QC, Enteritidis confirmation, core-SNP phylogenomics, recombination filtering, clade analysis, assembly/pangenome analysis, and AMR/virulence/plasmid analyses.

The complete original Part A notebook is preserved under provenance/ and is the authority for Part A statistical analyses. See analysis/part_a/README.md for the original-code runner and population separation.

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

### Part A — reproduce the complete ORIGINAL Rmd (33,831 full master)

**Authoritative source:** `provenance/SEN_Genomes2_uploaded_original_20261021.Rmd` (all 30 chunks). Do not substitute the modular 3,434-genome marker/source/diversity scripts for full-master analyses.

- Original Rmd chunks **3–22** run on the **33,831-genome master**, including phenotype trends, marker frequencies, International/Other-US source refinement, temporal graphs and diversity.
- Historical selection-generator chunk **23 is not rerun**: the original 3,434 selected genomes are a fixed final cohort.
- Original Rmd chunks **24–30** use that **3,434-genome cohort after QC merge**, while the master object remains 33,831 genomes.
- Deliberate edits only: strict CSV parsing and input paths; historical **≤2000** bin includes all dated pre-1990 genomes (33,508 dated total); consistent dot-plot source-facet Ns; Top-100 marker count axis to 25,000; consistent `SEN3` → `SEN` naming in the generated working copy; selection/QC audit.
- Full-master output: `results/part_a/full_master_original/`. Selected original output: `results/part_a/selected_original/`. QC merge output: `results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv`.

```bash
bash scripts/run_part_a.sh 2>&1 | tee original_part_a_rerun.log
```

**Additional October 2026 Part A2 outputs:** A Top-100 marker figure with axis extending to **25,000** but no endpoint tick label; paired full-master and selected-cohort Top-5 + Others + QS/FQS temporal count/percentage plots (exactly one category per genome); `SEN_33831_full_enriched_metadata.csv` (all original master columns and every original Rmd-derived field) and `SEN_3434_selected_enriched_QC_metadata.csv` (all original selected metadata plus QC and matched master fields). See `results/part_a/reconciled_metadata/` and `results/part_a/temporal_top5/` after the rerun; neither output replaces the archival data or original publication analyses.

See **[Updated Part A2 release and execution](docs/PART_A2_UPDATED.md)** and
**[Part A2 analysis details](analysis/part_a/README.md)**. This release tracks
only Part A2 source and cohort analyses; Table S5 verification is a separate
clade-analysis workstream. The legacy files `01_phenotype_trends.R`, `03_genotype_marker_summary.R`, `04_source_profiles.R`, and `05_diversity_composition.R` are *not* used by this official full-master runner because they had been executed on the wrong genomic denominator. If the reproduced International/Other-US source subgroups differ from an earlier manuscript Table S1, the validator reports **[REVIEW REQUIRED]** and saves the isolate-level audit without changing source assignments.

### Part B1 — read retrieval, trimming and quality control

The final Part A cohort provides the 3,434 selected accessions entering Part B.

1. Download Illumina reads with SRA Toolkit.
2. Trim and quality-filter paired-end reads with fastp.
3. Summarize raw and trimmed read quality with FastQC and MultiQC.
4. Classify reads with Kraken2 and retain Enterobacteriaceae reads (TaxID 543).
5. Require post-filter sequencing depth >=30x.

### Part B2 — species confirmation and typing

1. SeqSero2 k-mer serotyping.
2. Rapid SKESA assembly followed by SISTR.
3. Retain genomes confirmed as *Salmonella* Enteritidis by both workflows.
4. MLST typing.

The final confirmed study cohort contains 3,306 study isolates.

### Part B3 — core-genome SNP and phylogenetic analysis

1. Map reads to SEN reference P125109 (NC_011294.1) with Snippy.
2. Generate the core-genome alignment with snippy-core.
3. Annotate variants with SnpEff.
4. Mask recombination with Gubbins.
5. Infer the maximum-likelihood phylogeny with IQ-TREE.
6. Validate the final 3,306 study isolates plus P125109 against canonical clade metadata.
7. Perform clade-enriched and clade-defining SNP analyses.

### Part B4 — assembly, annotation and pangenome analysis

1. Generate de novo draft assemblies with Shovill/SPAdes.
2. Annotate assemblies with Prokka.
3. Reconstruct the pangenome with Panaroo in strict mode using MAFFT.
4. Test accessory-gene clade enrichment.
5. Define clade-associated accessory genes as >=90% within clade and <=5% outside clade.

### Part B5 — AMR, virulence and plasmid analysis

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

Part A validates original input contracts, then independently audits full 33,831-genome original figures/source counts and the 3,434-genome selected/QC handoff before reporting PASS. Part B includes a representative smoke
workflow and canonical metadata/clade-count checks.

## Provenance

The original working genomics repository is ajulojays/SENBio-Final.

The full, user-uploaded 30-chunk original is tracked at:

    provenance/SEN_Genomes2_uploaded_original_20261021.Rmd

The earlier 26-chunk `provenance/SEN_Genomes2.Rmd` remains historical provenance. The official Part A workflow extracts and executes the 30-chunk original with only documented deliberate edits, and preserves the frozen selected cohort rather than rerunning historical subsampling.

## Release status

This repository contains the complete manuscript-facing SEN workflow implementation spanning Part A and Part B.

A fresh run uses the workflow code, canonical data/metadata, and pinned software/database definitions stored in this release.
