# Publication archive and supplementary-data status (2026-11-02)

## Repository

Public URL: https://github.com/ShahShringiLab/SEN_genomics_pipeline
Default branch: main
This repository includes SEN Part A1 phenotype and Part A2 genome-selection/genotype analyses. The original historical R Markdown remains at `provenance/SEN_Genomes2_uploaded_original_20261021.Rmd`; the documented working notebook is `provenance/SEN_Genomes2.Rmd`.

**Caution:** The moving `main` branch is not a fixed archive. Cite a frozen Git commit/tag or a Zenodo DOI assigned to a release. No DOI should be stated until a Zenodo record exists.

## Supplementary raw data completeness — not yet confirmed

The repository presently carries input datasets and selection/cohort files, notably:
- `data/part_a/genotype/US-SEN_12192025.csv`
- `data/part_a/phenotype/FDASEN_AMRPhenotype2.csv`
- `data/part_a/selection/Finalupdated_SEN_Genomes1.csv` (archived selected cohort)
- `data/part_a/selection/Selected_US-SEN_12192025_QCAUDIT_R1_qr_clade.csv`
- `metadata/SEN_Genomes.csv`

The *publication supplementary raw data files* are separate submission artifacts. Their complete 33,831-genome export (original NCBI metadata + R-derived fields + audit + selected-QC annotations) and any other supplementary raw data attachments have **not been confirmed in this GitHub tree**. The 33,831-row full enriched CSV is generated locally by `analysis/part_a/09_temporal_top5_and_enriched_metadata.R`, which writes `results/part_a/reconciled_metadata/SEN_33831_full_enriched_metadata.csv`.

Before asserting completeness, confirm each submitted supplementary raw data file exists and opens; verify the expected full-master 33,831 unique genome rows, the submitted column layout and annotations, selected 3,434 pre-QC and 3,306 post-QC membership, absence of duplicated accessions, and table/figure linkage. Record checksums and verify contents of any repository archive or Zenodo upload. Missing QC/clade values for unselected genomes may be expected, not necessarily incomplete data.

## DOI archive procedure

1. Confirm the final analysis scripts and supplementary data files are correct.
2. Connect the public organization repository through https://zenodo.org/account/settings/github/ (organization owner may need to authorize Zenodo).
3. **After the connection is enabled**, publish a GitHub release on the final commit/tag.
4. Verify Zenodo ingested the release and use the exact DOI from Zenodo. A GitHub release alone does not provide a DOI.
5. If supplementary files are not included in the GitHub release, deposit them separately and cite their archived record, rather than claiming they are inside the code release.

## Data Availability Statement (before DOI issuance)

The analysis code and repository-hosted datasets are publicly accessible at https://github.com/ShahShringiLab/SEN_genomics_pipeline. A fixed revision is available through the repository's commit history. [Insert the final versioned release URL and DOI only after publication of the release and successful archival.] Supplementary raw data files accompanying the manuscript should be referenced individually and their completeness verified before submission.
