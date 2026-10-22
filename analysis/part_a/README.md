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
4. Original `FQR_marker_counts_top100` horizontal count-axis maximum changed from **22,000 to 250,000**, with the **250,000 tick label suppressed** (labeled ticks stop at 200,000). Marker values unchanged. This tenfold axis expansion is per the October 22 request and will visually compress the bars.

### Additional reconciled metadata and matched temporal stacks (October 22)

The **original Rmd and prior manuscripts' figures remain unmodified** by this additive deliverable. The complete original full-master analysis now exports every Rmd-created non-raw field (`SEN_33831_original_Rmd_derived.csv`) after original chunk 22. `analysis/part_a/09_temporal_top5_and_enriched_metadata.R` then creates:

- `results/part_a/reconciled_metadata/SEN_33831_full_enriched_metadata.csv`: **33,831 genomes**; preserves each of the 51 corrected original metadata columns, appends all `Rmd_` fields (primary/final source, year, marker variables, etc.), consistent `Audit_` columns and `QC_` status for selected genomes; QC entries remain blank for unselected genomes.
- `results/part_a/reconciled_metadata/SEN_3434_selected_enriched_QC_metadata.csv`: **3,434 selected genomes**; preserves all 69 original selection columns and all 9 appended QC columns (78 total), then joins exactly the same `Rmd_` and `Audit_` columns from its 33,831-genome master counterpart by `Run + Isolate`. The 3 corrected QR calls remain traceable through archived and corrected fields.
- `SEN_3434_master_selected_field_discrepancies.csv`: only discrepancies in source, collection year and QR status between original selected metadata and corresponding full-master records, with original values retained.
- `SEN_metadata_reconciliation_summary.csv`: counts and validations.

New charts in `results/part_a/temporal_top5/` render the **total master and selected cohorts separately**, each with count and percentage plots in PNG/TIFF/PDF. They use the same seven mutually exclusive categories: **QS/FQS, gyrA_D87Y, gyrA_D87N, gyrA_S83Y, gyrA_D87G, qnrB19 and Others**. Periods: **Pre-2018, 2018–2020, 2021–2025**. The master chart contains **33,508 dated genomes**; the selected chart contains **3,434**. Unlike plots of overlapping marker occurrences, **each genome is counted exactly once**, prioritizing the five markers in listed order when multiple are present. A dedicated `SEN_master_multimarker_bucket_priority_audit.csv` documents those rare assignments. This avoids totals exceeding cohort sizes.

No sampling is repeated, no original code is rewritten and no biological source result is forced to match the older Table S1. The independent validator `scripts/validate_enriched_metadata_outputs.py` verifies exact key and QC reconciliation and chart totals.

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

## Recover from an error after full-master analysis has completed

If original chunks 3–22 already ran and `results/part_a/full_master_original/SEN_33831_original_Rmd_derived.csv` exists, **you do not need to recalculate 33,831-genome full-master figures**. Reuse the verified original raw and Rmd-derived snapshots to rehydrate the complete `SEN` object and rerun only the selected-genome analyses (chunks 24–30), then the new metadata and paired temporal plots:

```bash
cd ~/SEN_PartA_clean
git pull --ff-only origin feature/part-a2-master-33831-parsed-snapshot-20261021
bash scripts/run_part_a_selected_resume.sh 2>&1 | tee selected_resume_rerun.log
```

The resume runner checks all prerequisite snapshots and exits rather than silently regenerating or substituting selected data. Its R log is `results/part_a/master_33831/original_execution/selected_resume.log`. It reconstructs the original binary `QR_genotype` alias omitted with historical selection chunk 23 and ensures `QR_status` is present in both plotting views, leaving the original full-master analysis unchanged. The older original marker-occurrence stacks may sum to more than genome totals when multiple marker hits occur; the separate seven-category `temporal_top5/` outputs count each genome once.

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
