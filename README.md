# Salmonella Enteritidis Genomics Pipeline

Reproducible publication workflow for the *Salmonella enterica* serovar Enteritidis (SEN) analyses described in the associated manuscript.

The repository implements the complete study workflow in two connected parts:

- **Part A** — phenotype resistance trends, QR/FQR genotype epidemiology, source/origin classification, marker diversity, and construction/validation of the final 3,434-genome cohort.
- **Part B** — read retrieval and QC, Enteritidis confirmation, core-SNP phylogenomics, recombination filtering, clade analysis, assembly/pangenome analysis, and AMR/virulence/plasmid analyses.

The original exploratory Part A notebook is preserved under provenance/. Linux-executable production scripts live under analysis/part_a/.

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

Input:

    data/part_a/phenotype/FDASEN_AMRPhenotype2.csv

The production analysis reproduces the Human-versus-Chicken temporal analysis of ciprofloxacin resistance, nalidixic-acid reduced susceptibility/resistance, and the combined QR/FQR phenotype metric with 95% confidence intervals.

Run:

    conda run --no-capture-output -n sen_stats \
      env SEN_ROOT="$PWD" \
      Rscript analysis/part_a/01_phenotype_trends.R

### Part A2 — genotype epidemiology and cohort construction

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
