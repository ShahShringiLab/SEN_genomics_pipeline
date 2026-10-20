# Part A data

This directory contains the committed data inputs and canonical cohort artifacts
used by the production Part A workflow.

## phenotype/FDASEN_AMRPhenotype2.csv

Aggregated Human and Chicken *Salmonella Enteritidis* quinolone /
fluoroquinolone phenotype counts by year. This is the input to
`analysis/part_a/01_phenotype_trends.R`.

Release contract:

- 58 rows;
- sources: Human and Chicken;
- one row per Year × Source;
- Ciprofloxacin-tested and DSC-tested denominators are identical in the
  committed dataset.

## genotype/US-SEN_12192025.csv

Upstream genotype/metadata table used for Part A2 source classification,
collection-year derivation, QR/FQR-marker analysis, and final cohort selection.

Release contract:

- 33,831 isolates;
- 51 columns;
- unique `Isolate` and `Run` values.

The source CSV is preserved unchanged. It contains known source-file quoting
artifacts: 24 literal escaped-comma sequences and one doubled-quote field around
`Exudate, Left Heel`. The production reader normalizes these only in memory
for parsing and then requires the expected dimensions.

## selection/Selected_updated_Isolate.csv

Final two-column (`Strain`, `Isolate`) selection key.

Release contract:

- 3,434 rows;
- 3,434 unique `Isolate` values;
- 3,434 unique `Strain` + `Isolate` key pairs.

## selection/Finalupdated_SEN_Genomes1.csv

Committed final Part A2 cohort and validation target.

Release contract:

- 3,434 isolates;
- unique `Isolate` and `Run` values;
- final source counts:
  - Human: 1,976
  - Chicken-US: 1,059
  - Other-US: 328
  - International: 71

Its Run accession set must match `metadata/SEN_Genomes.csv`, which is the
canonical Part B starting cohort.

## QR genotype audit note

Reapplying the notebook's original QR/FQR genotype rule to the upstream AMR text
yields three `gyrA_D87Y`-positive records whose committed final
`QR_genotype` is 0:

- PDT001102834.1
- PDT001168391.1
- PDT001168656.1

The production workflow preserves the archived values, writes the three-record
discrepancy audit to `results/part_a/cohort/QR_genotype_corrections.csv`, and
uses the curated AMRFinderPlus-core-derived QR status for production summaries:
2,288 positive and 1,146 negative.
