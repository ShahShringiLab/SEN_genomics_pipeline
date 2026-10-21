# Part A — phenotype and genotype/cohort analysis

Part A contains the Linux-executable implementation extracted from the original
SEN analysis notebook. The original notebook is retained separately as
provenance; the scripts in this directory are the reproducible publication
workflow.

## A1. Phenotype resistance trends

Script: analysis/part_a/01_phenotype_trends.R

Input:

    data/part_a/phenotype/FDASEN_AMRPhenotype2.csv

Outputs include the Human-versus-Chicken quinolone/fluoroquinolone phenotype
trend table, 95% confidence intervals, and publication PNG/PDF figures.

The script preserves the notebook rule excluding 2024 and 2025. It also checks
that Cip_Tested and DSC_Tested are identical before using the shared denominator
defined in the original analysis.

## A2. Genotype trends and final cohort construction

Script: analysis/part_a/02_build_final_cohort.R

Inputs:

    data/part_a/genotype/US-SEN_12192025.csv
    data/part_a/selection/Selected_updated_Isolate.csv
    data/part_a/selection/Finalupdated_SEN_Genomes1.csv

The script rebuilds the 3,434-isolate selected cohort from the master genotype
table and final two-column selection key. It reproduces collection-year,
QR/FQR-marker and source/origin derivations and verifies the rebuilt fields
against the committed final 3,434-genome table.

Script: analysis/part_a/03_genotype_marker_summary.R

Produces marker-frequency, marker-burden, multi-marker-combination,
source-by-year QR/FQR and top-marker prevalence outputs. Marker-frequency and
composition outputs are token-based, following the notebook marker definition;
binary QR/FQR status summaries use the committed final `QR_genotype` after the
three-record discrepancy audit.

Script: analysis/part_a/04_source_profiles.R

Produces final source counts, food-origin summaries, refined Other-source
animal categories, and International/Other-US animal-source × QR/FQR tables.

Script: analysis/part_a/05_diversity_composition.R

Produces Hill-number diversity estimates (q0, q1, q2), bootstrap 95% confidence
intervals, QR/FQR composition, and global/pairwise permutation tests.

The original notebook used manually supplied plot letters for one diversity
figure. The production script does not hard-code significance letters; it writes
the permutation-test tables directly.

## Run

From the repository root:

    bash scripts/run_part_a.sh

The runner is fail-fast. It validates canonical inputs before analysis, executes
the five modular R stages in order, and then runs
`scripts/validate_part_a_outputs.py`. It prints `SEN PART A: PASS` only when
all scientific/data contracts and required output artifacts pass.

All generated Part A outputs are written beneath results/part_a/.


## Release safeguards

The master genotype CSV is kept unchanged in `data/part_a/genotype/`. The
production reader performs two parser-only normalizations in memory for known
source-file CSV artifacts, then requires exactly 33,831 rows and 51 columns.

Production QR/FQR status is derived from the AMRFinderPlus `AMR genotypes core`
field. This field provides the stable curated determinant calls and avoids a
known source-CSV quoting artifact that shifted the broader AMR field in three
records: `PDT001102834.1`, `PDT001168391.1`, and `PDT001168656.1`.
All three contain `gyrA_D87Y=POINT` in the core field and are therefore
corrected to QR/FQR marker-positive in production. The archived values are kept
as `QR_genotype_archived` for provenance.

The rebuilt Part A cohort also joins `metadata/final_clade_metadata.tsv` by
SRA `Run`. Exactly 3,306 genomes receive final clade labels; 128 selected
genomes were excluded before the final phylogeny and retain a missing `Clade`
value with `Final_phylogeny_status = "Excluded before final phylogeny"`.

A compact committed audit table containing Run, Isolate, original/corrected QR
status, detected QR markers, clade, phylogeny inclusion status, and QC audit is:

    data/part_a/selection/Selected_US-SEN_12192025_QCAUDIT_R1_qr_clade.csv

No script installs packages at runtime. Dependencies are supplied by the pinned
`sen_stats` Conda environment. The archived Rmd is provenance only; production
execution does not depend on notebook chunk state, interactive R sessions, or
workstation-specific paths.


## October 2026 review note

The lightweight 3,434-row Run/Isolate/Clade sidecar matches the curated
3,434-row QR/clade audit (including 128 outside the final phylogeny). Three
`gyrA_D87Y`-positive Clade 3 isolates have QS/FQS curated overrides in the
newer clade-summary code. These are different analytical constructs and
must not be conflated with absence of a QR/FQR molecular determinant.
