# SEN Part A — reproducible original analysis (October 2026)

## Analytical authority and population boundaries

**The complete 30-chunk uploaded original Rmd is the scientific source of truth:**
provenance/SEN_Genomes2_uploaded_original_20261021.Rmd.

The main command is `bash scripts/run_part_a.sh`.

- **Original chunks 3–22:** Phenotypic trends, marker frequencies and burden, source refinements, International/Other-US Table S1, temporal dot and stacked plots, and Hill-number diversity. The genomic input is always **the full 33,831-genome SEN master**, never the 3,434-genome selected metadata. Phenotype data are read from their original surveillance file.
- **Original chunk 23:** Historical, stochastic selection generation **is not rerun**. The original final 3,434-genome selection has already been frozen and committed.
- **Original chunks 24–30:** Downstream selected-genome analyses use the fixed **3,434-genome QC-merged metadata**; `SEN` remains the independent 33,831-genome master object. Input file paths are updated but analyses are preserved.
- **Independent QC reporting:** The validated 3,434 selected metadata, 3,306 phylogeny genomes and 128 exclusions, three molecular-QR corrections, 2018–2025 selection and QC tables are additional audits, not replacements for original full-master analyses.

### Only deliberate modifications to the original Rmd

The original committed file remains unchanged; `scripts/prepare_original_part_a2.py` generates a working copy with:
1. Input paths redirected to the strictly parsed 33,831 × 51 master snapshot and committed phenotype CSV.
2. Four/six unintended `Collection_year >= 1990` exclusions removed in temporal analyses. All valid years **≤2000** remain in the combined historical bin, so **33,508** dated genomes contribute, and **323** undated genomes are excluded.
3. Original temporal facet **N** restricted to the same eligible dated genomes as the marker-percentage denominators. No marker percentages or methods are changed by this label fix.
4. Original `FQR_marker_counts_top100` horizontal count-axis maximum changed from **22,000 to 25,000**. Marker values and styling unchanged.

### Source decomposition audit

The original master's manuscript table `SEN_TableS1_final_source_decomposition.csv` must contain:

| Group | Subsource | Genome n |
|---|---|---:|
| International | Others | 67 |
| International | Chicken | 4 |
| Other-US | Others | 448 |
| Other-US | Cattle | 60 |
| Other-US | Turkey | 48 |
| Other-US | Swine | 21 |
| Other-US | Sheep/goat | 5 |

Expected group totals: **International 71**, **Other-US 582**. Corrected parser changes should be **audited, never forced to preserve old numbers**. First rerun of the complete original code found `Other-US: Swine=43, Others=426` (versus manuscript `Swine=21, Others=448`). Cattle=60, Turkey=48 and Sheep/goat=5 are unchanged. Source isolate-level classification will be audited; do not silently treat 43 as an updated biological truth until its 22-record discrepancy is inspected. Original Rmd chunk 18 also uses two obsolete group labels; the generated execution script redirects its filter to the already-derived `Final_source` to avoid empty charts, leaving the original animal-source and QR rules unchanged. The end-of-run validator prints **[REVIEW REQUIRED]** and records all differences in `REPRO_AUDIT.json` rather than silently changing Table S1; it does not falsely fail otherwise valid computations over a historical manuscript reference. A 653-isolate source audit and an isolate-level Swine audit are exported for reconciliation.

### Outputs (never mix denominators)

- `results/part_a/full_master_original/`: The **full 33,831** original marker/source/diversity/temporal tables and figures, including original Table S1, `SEN_fullmaster_source_isolate_audit.csv`, `SEN_OtherUS_swine_isolate_audit.csv`, and a machine-readable reproducibility audit.
- `results/part_a/selected_original/`: Original selected-genome figures and tables from Rmd chunks 24–30.
- `results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv`: Complete **3,434 × 78** metadata with all original 69 fields and nine QC/provenance fields.
- `results/part_a/master_33831/00_parsed_raw/`: Strict 33,831 × 51 snapshot, field-alignment checks and checksums.
- `results/part_a/master_33831/sampling_qc_table/`: Master-to-selection sampling weights/QC table.
- `results/part_a/master_33831/original_execution/run.log`: Full original execution log.

The alternative extracted scripts `01_phenotype_trends.R`, `03_genotype_marker_summary.R`, `04_source_profiles.R`, and `05_diversity_composition.R` are **not invoked by the publication reproducer**, because their 3,434-genome inputs would mislabel full-master results. Their old outputs must not be used for master population claims.

## Fresh clone and rerun

Keep the previous working folder until the new clone passes.

```bash
git clone --single-branch --branch feature/part-a2-master-33831-parsed-snapshot-20261021 \
  https://github.com/ShahShringiLab/SEN_genomics_pipeline.git \
  ~/SEN_PartA_clean
cd ~/SEN_PartA_clean
conda env list | grep sen_stats
bash scripts/run_part_a.sh 2>&1 | tee original_part_a_rerun.log
```

If the Conda environment is not already installed, create it using the repository environment files before rerunning.

**Do not delete previous files until all PASS checks finish.** Because the 33,831-genome raw CSV and original Rmd are committed, the fresh clone should not require manually copying inputs. Runtime analyses have to execute on a machine with the `sen_stats` environment and R packages; the GitHub commit alone does not establish successful execution.

### Original Rmd execution compatibility (not new biological methods)

The uploaded notebook assumes an interactive R session. The runner reconstructs missing `SEN_plot` and `SEL_plot` views for chunk 25 using the original three periods and original binary marker variable; it never resamples the fixed cohort. Chunk 19's `src_N` deliberately discards its intermediate `N_total` after generating source labels, so the validator checks `den_year`, the eligible full-master year count and the retained `Source_lab` rather than referencing a nonexistent column. Chunk 18's obsolete source-group filter is translated to the existing `Final_source` groups `International` and `Other-US` in the execution copy only, to prevent empty results.
