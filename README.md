# Salmonella Enteritidis Genomics Pipeline

Reproducible analysis workflow for the *Salmonella enterica* serovar Enteritidis (SEN) genomic analyses described in the associated manuscript.

This repository is organized as a publication-facing pipeline rather than a historical working directory. It separates scientific analysis scripts, software environments, configuration, reference documentation, and validation tests.

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

Scientific thresholds are centralized in `config/pipeline_config.yaml`. Machine-specific paths should not be committed. Use `config/paths.example.env` as the template for local/HPC paths.

## Validation

The first validation target is a smoke test using one representative genome from each major SEN clade plus P125109. See `tests/smoke_test/README.md`.

## Provenance

The original working analysis repository is `ajulojays/SENBio-Final`. This repository is the cleaned, publication-facing implementation.
