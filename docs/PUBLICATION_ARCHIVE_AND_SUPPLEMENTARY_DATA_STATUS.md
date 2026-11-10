# SEN publication source and data availability

## Canonical repository and release

Public repository: https://github.com/ShahShringiLab/SEN_genomics_pipeline

The unified publication workflow is organized under `workflow/A1_Phenotype/`, `workflow/A2_Genotype/` and `workflow/B1_Reads_QC/` through `workflow/B5_AMR_Virulence_Plasmids/`, with filenames retained. The README links each scientific script; `docs/FULL_RUN_FRESH_LINUX.md` records the complete execution sequence.

The intended publication source snapshot is **v1.0.0**, with a GitHub Release containing a tagged `.tar.gz` source archive and a SHA-256 checksum file. The definitive citation should include the tag and exact Git commit SHA. A tag and GitHub Release provide a fixed version in the repository, but are not a separate DOI-backed preservation service. **Do not claim a DOI.**

Release page (available only after the release workflow succeeds):
https://github.com/ShahShringiLab/SEN_genomics_pipeline/releases/tag/v1.0.0

## Included canonical data and analysis provenance

- The phenotype and genotype/source inputs are in `data/part_a/`; the full master contains 33,831 genomes.
- The fixed selection contains 3,434 genomes, with 3,306 included after QC and 128 excluded.
- The original exploratory R Markdown is retained under `provenance/` for method provenance.
- The final clade metadata are in `metadata/final_clade_metadata.tsv`, with 3,306 study isolates plus one P125109 reference.
- Historical Other-US swine count **43** was independently confirmed during publication review; any earlier manuscript number of 21 is superseded. Source counts must remain traceable to the source audit rather than silently rewritten.
- Raw FASTQs, larger intermediate outputs, the Kraken2 database and other external scientific assets are not part of the Git source archive. Their retrieval and pinned versions are documented in `config/` and `docs/FULL_RUN_FRESH_LINUX.md`.

## Validation and interpretation

Source-contract checks, the Part A enriched metadata checks, and representative Part B execution were performed during migration. Contract validation is not equivalent to independent confirmation of every historical manuscript-scale result. Software-version and external-service provenance qualifications are retained in `docs/publication_audit.md`.

## Data availability wording

The genomic metadata analyzed in this study were obtained from the NCBI Pathogen Detection Isolates Browser. The complete analysis scripts, canonical repository-hosted inputs, and reproducibility documentation are available at https://github.com/ShahShringiLab/SEN_genomics_pipeline. The fixed publication software version is identifiable by the `v1.0.0` release tag and its Git commit identifier. Supplementary raw data files accompany the manuscript. No independent code DOI is claimed.
