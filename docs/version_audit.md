# Software version audit

This document separates **current workstation observations** from the
**manuscript Methods specification**. They are not assumed to be identical.

## Current workstation observations used to seed environment files

| Stage | Environment | Observed package/version |
|---|---|---|
| Reads/QC | senbio | fastp 1.0.1; FastQC 0.12.1; MultiQC 1.14; SeqKit 2.12.0 |
| SRA | external / seqsero2_env | SRA Toolkit 3.2.1 observed in seqsero2_env |
| Kraken | kraken2_env | Kraken2 2.1.3 |
| Serotyping | seqsero2_env / sistr_env | SeqSero2 1.3.1; SKESA 2.5.1; SISTR 1.1.3 |
| MLST | mlst_env | mlst 2.33.1 |
| Core SNP | snippy_env | Snippy 4.6.0; SnpEff 5.1 |
| Recombination | gubbins_env | Gubbins 3.4.3; VeryFastTree 4.0.5 |
| Phylogeny | iqtree_env | IQ-TREE 2.2.6 |
| Assembly | shovill_env | Shovill 1.0.9; SPAdes 4.0.0 |
| Annotation/pangenome | prokka_env / panaroo_env | Prokka 1.14.6; Panaroo 1.6.0 |
| AMR | amrfinder_env | NCBI AMRFinderPlus 4.2.5 |
| ABRicate | abricate_env | ABRicate 1.2.0 |
| R | r_env | R 4.5.2 |

## Manuscript-facing discrepancies requiring resolution

The current workstation state is not automatically the historical execution
state. In particular:

- Methods specify Kraken2 2.1.3, which matches `kraken2_env`; the separate
  `contamination_qc` environment currently contains a newer Kraken2.
- Methods specify MLST 2.23.0, while the current `mlst_env` contains 2.33.1.
- Methods specify Gubbins 3.0, while the current `gubbins_env` contains 3.4.3.
- Methods specify IQ-TREE 2.2.6, which matches `iqtree_env`; the IQ-TREE copy
  bundled in `gubbins_env` is newer and should not be used for the final-tree
  reconstruction.
- Methods describe AMRFinderPlus v3, while the current `amrfinder_env` contains
  4.2.5.

Before publication, each mismatch must be resolved by either recovering the
historically executed version from provenance/logs or updating the Methods if
the analysis was in fact run with the currently observed version.

The environment YAMLs in this repository are therefore **reconstruction
starting points based on observed workstation packages**, not proof of the
historical software versions used for every final result.
