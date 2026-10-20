# Part A on a fresh Linux machine

Part A covers the epidemiologic phenotype analysis and the genotype/source/cohort
analysis that precede the comparative-genomics workflow.

## 1. Clone the repository

    git clone https://github.com/ShahShringiLab/SEN_genomics_pipeline.git
    cd SEN_genomics_pipeline

For a published analysis, checkout the corresponding release tag or commit.

## 2. Create the declared environments

    bash scripts/setup_all_envs.sh

Part A uses the sen_stats environment. The environment includes R, tidyverse,
lubridate and patchwork.

## 3. Verify Part A inputs

The repository must contain:

    data/part_a/phenotype/FDASEN_AMRPhenotype2.csv
    data/part_a/genotype/US-SEN_12192025.csv
    data/part_a/selection/Selected_updated_Isolate.csv
    data/part_a/selection/Finalupdated_SEN_Genomes1.csv

The master genotype file is the upstream metadata/genotype table. The
Selected_updated_Isolate.csv file is the final 3,434-isolate selection key.
Finalupdated_SEN_Genomes1.csv is the committed final Part A2 cohort used as
the validation target and handoff to Part B.

## 4. Run Part A

    bash scripts/run_part_a.sh

A successful run ends with:

    SEN PART A: PASS

The runner executes, in order:

1. strict Part A input and Part A -> Part B handoff validation;
2. phenotype trend analysis;
3. final 3,434-isolate cohort reconstruction and validation;
4. QR/FQR marker summaries and temporal/source trends;
5. source/origin and International/Other-US animal-source profiling;
6. QR/FQR diversity and composition analysis;
7. generated-output validation.

The final validator checks the rebuilt cohort size/uniqueness, all cohort
reconstruction mismatch counts, expected final source counts, phenotype result
dimensions, diversity source coverage, and the presence of required figure
artifacts. Any failed check stops the workflow before the PASS banner.

Generated files are written beneath:

    results/part_a/

The original `provenance/SEN_Genomes2.Rmd` is not executed during production
runs. It is retained as the analytical provenance record from which the modular
Linux scripts were extracted.

## 5. Scientific contracts

Part A1 preserves the original phenotype analysis contract, including the
2024/2025 exclusion and normal-approximation 95% confidence intervals. The
script explicitly verifies the shared Ciprofloxacin/DSC denominator before
computing the reported percentages.

Part A2 preserves the notebook's broad QR/FQR marker definition, source
classification logic, collection-year extraction, FDA/CFSAN food-origin logic,
and final 3,434-isolate selection key. The archived master CSV is never edited;
known source-file quoting/escaped-comma artifacts are normalized only in memory
for parsing. The reconstructed cohort is compared against the committed final
cohort before downstream Part A2 analyses continue.

Three upstream gyrA_D87Y-positive isolates differ from the committed final
QR_genotype. They are retained in an explicit discrepancy audit. Production QR/FQR binary calls
are instead recomputed from the curated AMRFinderPlus `AMR genotypes core` field,
with the archived values preserved for provenance.

## 6. Handoff to Part B

The Part A2 final cohort contains the 3,434 selected genomes from which Part B
performs sequence retrieval, QC, Enteritidis confirmation and subsequent
comparative-genomics analyses. Part B ultimately retains the final 3,306 study
SRRs plus P125109 for downstream phylogenomic analyses.
