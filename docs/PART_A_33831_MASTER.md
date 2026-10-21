# Part A2 master-population run (33,831 genomes)

The original Part A2 notebook remains unchanged at `provenance/SEN_Genomes2.Rmd`. **The 33,831-genome master population is the denominator** for original-population QR/FQR, source, marker and temporal summaries. The 3,434 selected and 3,306 final-phylogeny cohorts are separate sampling/validation stages.

## Clone and run on Linux / WSL Ubuntu

```bash
source ~/miniconda3/etc/profile.d/conda.sh
cd ~
git clone https://github.com/ShahShringiLab/SEN_genomics_pipeline.git SEN_master_A2
cd SEN_master_A2

# First create sen_stats if needed; see environments/12_stats.yaml.
conda env create -f environments/12_stats.yaml
conda install -n sen_stats -c conda-forge r-writexl -y

bash scripts/run_part_a_master.sh
```

If `sen_stats` already exists, **skip the environment creation** and install `r-writexl` only if needed. A previously cloned checkout can use `git pull origin main`, but preserve local modifications before pulling.

## Workflow and provenance

1. `analysis/part_a/00_parse_master_33831.R` reads `data/part_a/genotype/US-SEN_12192025.csv` through the original `read_sen_master_csv()`. That function normalizes 24 literal escaped commas and one malformed quoted field **in memory**. It strictly checks 33,831 records × 51 original columns, unique Isolate/Run, and a basic BioProject/BioSample/Run field-alignment audit.
2. Before any new variables are added, it saves `results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.xlsx` and `.csv`, plus alignment audit and file checksums. The original tracked CSV remains unchanged.
3. `scripts/run_part_a.sh` executes original A1 FDA phenotype trends and the 3,434-isolate cohort A2 analysis, including input/output validators. `scripts/validate_part_a_outputs.py` correctly recognizes the 128 missing/NA clade labels as excluded.
4. `analysis/part_a/06_master_33831_analysis.R` reads **only the saved parsed 51-column CSV** before adding year, source, origin and QR-marker fields. The QR call comes from `AMR genotypes core`. It audits selection (3,434), final phylogeny (3,306), missing/out-of-window years, and (when available) compares the 3,434 overlapping rows with the rebuilt selected cohort.

## Output locations

- **Original 51-column spreadsheet:** `results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.xlsx`
- **Derived master spreadsheet:** `results/part_a/master_33831/SEN_33831_DERIVED_ANALYSIS_READY.xlsx`
- **Derived master CSV:** `results/part_a/master_33831/SEN_33831_DERIVED_ANALYSIS_READY.csv`
- **Population-source and year eligibility:** `master_source_temporal_qc_audit.csv`
- **Marker frequency:** `master_marker_frequency.csv`
- **Marker burden:** `master_marker_burden_per_isolate.csv`
- **Annual QR/FQR proportions:** `master_QRFQR_by_source_and_year.csv`
- **Year/source/QR sampling weights:** `master_selection_weights_year_source_marker.csv`
- **Figures (PNG + PDF):** annual and five-year marker dot plots; stacked QR/FQR percentages; stacked counts; available/selected/analyzed trends; sampling coverage.

For master temporal charts, dates outside **1990–2025** and missing dates are excluded from time-specific N. Years up to 2000 are grouped into `≤2000` on grouped figures; within-year percent denominators are computed separately by source/year. The master population size remains 33,831.

## Scope and limitations

This implements the source-population parsing and the principal population-level temporal/marker/coverage figures, **not all exploratory chunks** of the 6,469-line original Rmd. In particular, notebook-specific refinements and full master-scale diversity, permutation and multivariate outputs have not yet been ported; do not represent this as a complete recreation of every legacy figure.

Do not modify the archived 51-column master data to repair the three shifted rows in the **historical selected file**. Those source data are normalized while reading; the corrected QR/FQR status is always recomputed from the AMRFinderPlus core genotypes, and differences with the selected cohort are audited.

The scripts and data were uploaded to GitHub; complete execution on the user's WSL machine is still required to validate new master outputs.
