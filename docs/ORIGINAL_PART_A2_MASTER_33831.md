# Original Part A2 — 33,831-genome master reproduction

The original population-level analysis is documented in `provenance/SEN_Genomes2.Rmd`. Do **not** substitute the selected 3,434 genomes as the denominator for 33,831-genome epidemiology.

## Clone and run

```bash
git clone --branch feature/part-a2-master-33831-parsed-snapshot-20261021 \
  https://github.com/ShahShringiLab/SEN_genomics_pipeline.git SEN_MasterA2
cd SEN_MasterA2
source ~/miniconda3/etc/profile.d/conda.sh
# Reuse your existing sen_stats environment; writexl is needed for Excel export.
conda install -n sen_stats -c conda-forge r-writexl -y
bash scripts/run_master_part_a2.sh
```

### Three stages
1. `analysis/part_a/00_parse_master_33831.R`: use `read_sen_master_csv()`, audit BioProject/BioSample/Run field alignment, then export the **unaugmented** 33,831-by-51 snapshot before adding variables. Does not modify the original CSV.
2. `scripts/run_part_a.sh`: run validated A1 phenotype plus A2 selected-cohort analyses and compare 3,434 selected versus 3,306 final-phylogeny genomes.
3. `analysis/part_a/06_master_33831_analysis.R`: derive collection year, QR/FQR status and source/origin using the canonical functions; cross-check selected and post-QC memberships; produce master population analyses and plots.

### Raw snapshot
`results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.xlsx`

Also: `SEN_33831_PARSED_51COL.csv`, `parsed_master_field_alignment_audit.csv`, `file_checksums.csv`.

### Derived outputs
Under `results/part_a/master_33831/`: analysis-ready CSV/XLSX, source/year QC audit, annual QR/FQR tables, marker frequencies/burden, two dot plots (year, five-year bins), stacked count and percentage plots, cohort-coverage plots, sampling weights, and master-vs-rebuilt discrepancy audit.

**Publication caveat:** This modular implementation restores key original Part A2 master-population analyses but does not assert that every exploratory chunk of the archived Rmd has been identically reproduced. The pipeline must run and audits must pass locally before publication.
