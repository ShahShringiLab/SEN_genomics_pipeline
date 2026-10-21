# Salmonella Enteritidis Genomics Pipeline

Reproducible analysis workflow for the *Salmonella enterica* serovar Enteritidis (SEN) genomic analyses described in the associated manuscript.

This repository contains the complete publication-facing implementation of the SEN genomic analysis workflow used for manuscript reproduction. It includes the full analysis code, pinned software environments, database/reference bootstrap utilities, canonical metadata, a four-genome end-to-end smoke test, downstream clade/SNP analyses, and fresh-Linux execution guides.

## Repository structure

```text
SEN_genomics_pipeline/
├── README.md
├── CITATION.cff
├── config/
├── environments/
├── workflow/
├── docs/
├── reference/
├── metadata/
└── tests/
```

## Reproduce the workflow on a fresh Linux machine

Two execution paths are documented separately.

### Four-genome smoke test

Use this first on a new machine. It validates the complete local workflow,
including downstream `analysis/clades` and `analysis/snps` runtime behavior.

See: `docs/SMOKE_TEST_FRESH_LINUX.md`

After installing Git and Conda, the principal command is:

```bash
bash tests/smoke_test/run_smoke_pipeline.sh
```

A successful run ends with:

```text
SEN END-TO-END SMOKE PIPELINE: PASS
```

### Full manuscript reproduction

The full workflow begins with the canonical 3,434-SRR metadata cohort and uses
the frozen historical final cohort of 3,306 study SRRs plus P125109 for the
publication-facing downstream analyses.

See: `docs/FULL_RUN_FRESH_LINUX.md`

The full guide contains copy-paste commands for environment creation, SRA
download, QC/typing, Kraken2, Snippy/SnpEff, Gubbins, IQ-TREE, Shovill,
Prokka, Panaroo, AMRFinderPlus/ABRicate, clade statistics, SNP analyses, and
the final publication audit. External EggNOG/PHASTER/HMMER steps are identified
explicitly rather than silently substituted.

## Analysis workflow

### 1. Read retrieval, trimming, and quality control

1. Download Illumina reads with SRA Toolkit.
2. Trim and quality-filter paired-end reads with fastp.
3. Summarize raw and trimmed read quality with FastQC and MultiQC.
4. Classify reads with Kraken2 and retain Enterobacteriaceae reads (TaxID 543).
5. Exclude genomes with post-filter sequencing depth <30×.

### 2. Species confirmation and typing

1. SeqSero2 k-mer serotyping.
2. Rapid SKESA assembly followed by SISTR.
3. Retain genomes confirmed as *Salmonella* Enteritidis by both workflows.
4. MLST typing.

### 3. Core-genome SNP and phylogenetic analysis

1. Map quality-filtered reads to SEN reference P125109 (NC_011294.1) with Snippy.
2. Generate the core-genome SNP alignment with snippy-core.
3. Annotate variants with SnpEff.
4. Mask recombination with Gubbins.
5. Infer the maximum-likelihood phylogeny with IQ-TREE.
6. Perform clade-enriched and clade-defining SNP analyses.

### 4. Assembly, annotation, and pangenome analysis

1. Generate de novo draft assemblies with Shovill/SPAdes.
2. Annotate Shovill assemblies with Prokka.
3. Reconstruct the pangenome with Panaroo in strict mode using MAFFT.
4. Test accessory-gene clade enrichment.
5. Define clade-associated accessory genes as ≥90% within clade and ≤5% outside clade.

### 5. AMR, virulence, and plasmid analysis

1. AMRFinderPlus with `--plus --organism Salmonella`.
2. ResFinder, VFDB, and PlasmidFinder through ABRicate.
3. Binary presence calls require ≥80% coverage and ≥90% identity.
4. Clade-stratified statistical analyses use omnibus tests, sparse-count handling, pairwise Fisher tests, and multiple-testing correction.

## Reproducible environments

Software environments are stored in `environments/`. Exact versions stated in the manuscript are pinned. Tools whose versions were not specified in the manuscript are intentionally left unpinned until the final lab-computer environment is reconciled.

## Configuration

Publication-facing scientific thresholds are encoded explicitly in the relevant workflow/analysis scripts and audited by `scripts/validate_analysis_contracts.py`. Frozen external database releases are recorded in `config/database_sources.env`. Machine-specific paths should not be committed; use `config/paths.example.env` as the template for local/HPC paths.

## Validation

The representative four-isolate smoke suite is documented in `tests/smoke_test/README.md` and `docs/SMOKE_TEST_FRESH_LINUX.md`. Full manuscript-scale execution is documented in `docs/FULL_RUN_FRESH_LINUX.md`.

## Provenance

The original working analysis repository is `ajulojays/SENBio-Final`. This repository is the cleaned, publication-facing implementation.


## Release status

This repository contains the complete manuscript-facing SEN pipeline implementation.

A fresh run from this repository uses the complete SEN workflow and the pinned software and database versions defined for the current release.

Run:

```bash
bash scripts/publication_audit.sh
```

to verify the repository's publication contracts on a fresh clone.
