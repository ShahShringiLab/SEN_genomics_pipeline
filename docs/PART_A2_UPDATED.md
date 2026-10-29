# SEN Part A2 — consolidated updated workflow (23 October 2026)

**Status:** Part A2-only update for review. This is the canonical execution and
documentation entrypoint for the revised Part A2 on this branch; **no Table S5
AMRFinder-by-clade statistical verification or manuscript Table S5 replication
belongs to this update**. Those analyses remain in the independent clade-analysis
workstream and are not called by the Part A2 runners.

## Scientific authority and population separation

The complete **30-chunk** uploaded original is retained, SHA-256 protected, at
`provenance/SEN_Genomes2_uploaded_original_20261021.Rmd`.
It is the original scientific reference, not the older 26-chunk archived notebook.

| Stage | Source / genome population | What it does |
|---|---|---|
| Original chunks 3–22 | **33,831** full-master SEN records | Original phenotype, QR/FQR genotype and marker, source/origin, temporal, diversity analyses and figures |
| Original chunk 23 | **Not rerun** | Historical stochastic selection; the validated **3,434** selection is frozen |
| Original chunks 24–30 | **3,434** selected records (after QC merge) | Original selected-cohort analyses, figures and tables |
| Part B handoff (informational) | **3,306 included / 128 excluded** | QC reconciliation; not a different Part A2 denominator |

Do **not** substitute the 3,434-genome cohort for full-master Part A2
analyses or rerun the random selection.

## What is updated in Part A2

The full analysis is generated from the original Rmd, with specific documented
working-copy compatibility changes:

1. **Input parsing and provenance:** `analysis/part_a/00_parse_master_33831.R`
   validates the full input as 33,831 records × 51 original fields, preserving
   original records and writing the parsed CSV and XLSX snapshot. The original
   uploaded 30-chunk source is never overwritten. `scripts/prepare_original_part_a2.py`
   redirects notebook inputs to that snapshot and the tracked FDA phenotype data.
2. **Temporal eligibility:** remove unintended pre-1990 cutoffs where the
   original chart already defines a ≤2000 bin. This retains **33,508 dated**
   master genomes for temporal analyses; **323 undated** are not used in dated
   graphs. Facet N is aligned with the plotted eligible population.
3. **Top-100 plot:** the `FQR_marker_counts_top100` axis extends to 25,000,
   with labels stopping at 20,000. Counts and scientific statistics are unchanged.
4. **Object consistency:** all five original `SEN3` references in the source-
   breakdown section of the **generated working copy** are changed to `SEN`
   (one assignment; four `plot_iso_counts` calls). The archived 26-chunk
   editable copy is also normalized. The SHA-locked 30-chunk upload remains
   unchanged. No source reclassification is introduced by this naming change.
5. **Original Rmd execution compatibility:** `scripts/build_original_a2_master_execution.py`
   extracts the original analysis chunks with 33,831/3,434 population guards,
   corrects an obsolete source-group label in its generated working script, and
   avoids rerunning historical chunk 23.
6. **Selected-cohort QC and resume:** `scripts/merge_qc_full_metadata.py`
   retains all **3,434** selected records, original fields and the QC audit.
   `scripts/run_part_a_selected_resume.sh` resumes only from validated
   master and selected snapshots, restoring needed original plotting aliases
   without resampling. Three QR corrections remain individually auditable.
7. **Additive epidemiological outputs:** `analysis/part_a/09_temporal_top5_and_enriched_metadata.R`
   generates separately scaled master/selected temporal count and percentage
   charts for Pre-2018, 2018–2020, and 2021–2025, with seven exclusive classes:
   QS/FQS, gyrA_D87Y, gyrA_D87N, gyrA_S83Y, gyrA_D87G, qnrB19, Others. It
   also writes full and selected enriched metadata with source and QC audit
   fields, preserving the original fields.
8. **International source breakdown:** `analysis/part_a/10_international_animal_source_plot.R`
   produces the independent International **N=71** decomposed-source plots
   and isolate-level Excel audit with original isolation-source wording retained.
9. **Within-scope tables:** Original Part A2 source-decomposition
   **Table S1**, and post-2018 QR/FQR marker-genotype **Table SS1** are retained
   in their original/extended Part A2 workflow and documentation. These are
   distinct from the **excluded Table S5** clade AMRFinder analysis.

## Canonical execution from a fresh or existing checkout

```bash
cd ~/SEN_PartA_clean
git fetch origin
git switch feature/part-a2-updated-20261023 || \
  git switch -c feature/part-a2-updated-20261023 \
    --track origin/feature/part-a2-updated-20261023

# Fast static check; checks source, SHA, extraction and cohort boundaries.
python3 scripts/check_original_rmd_contract.py

# Full original Part A2 and selected follow-ups (no Table S5 operations).
bash scripts/run_part_a.sh 2>&1 | tee updated_part_a2_run.log
```

`sen_stats` must be installed with the repository-declared R and Python
packages (see `environments/12_stats.yaml`). The full 33,831-record R run
must be performed locally; GitHub's fast checks validate structure and
provenance, **not the entire scientific computation**.

If full master chunks 3–22 have already completed successfully and the
validated snapshots exist, avoid an expensive rerun:

```bash
bash scripts/run_part_a_selected_resume.sh 2>&1 | tee part_a2_selected_resume.log
```

## Final-run CSV parsing and QR change audit

Every full run now writes `results/part_a/master_33831/00_parsed_raw/SEN_CSV_Malformed_Before_After.xlsx` with exactly two sheets:

- **Before:** the affected source CSV rows as parsed *without* the in-memory quoting/escaped-comma fixes, including raw field counts, row positions and displaced values.
- **After:** the corresponding records from the validated 33,831 × 51 corrected master snapshot, in the original column layout.

Additional audit files in the same folder:

- `SEN_CSV_Parsing_Field_Differences.csv` — affected CSV row, before/after field values, and corrected Run/Isolate identifiers. The 24 escaped commas and one doubled-quote anomaly are tracked as source occurrences; they are **not** assumed to be 25 distinct isolates.
- `SEN_CSV_Parsing_Change_Summary.csv` — counts of affected rows, changed cells and source-text anomalies.
- `SEN_QR_Genotype_Before_After_3_Genomes.csv` — each of the **three** audited changes between original and corrected selected-cohort QR genotypes, retaining accession IDs and marker information.

The 22-genome Swine/Others difference was traced by manual isolate review to the CSV parsing issue. The source audit keeps the raw and corrected values so this attribution can be checked at record level; this does **not** mean the 22 genomes changed full-master source-group totals. Corrected QR calls may affect QR-based counts and plots, but **no numerical effect is asserted without comparing outputs**. The original temporal-plot population remains **33,508 genomes with collection years** (323 missing-year records excluded), and this audit does not modify any figure or statistical algorithm.

## Main outputs

| Directory | Deliverable |
|---|---|
| `results/part_a/master_33831/00_parsed_raw/` | Strictly parsed 33,831 × 51 raw snapshot, XLSX/CSV, alignment audit and checksums |
| `results/part_a/full_master_original/` | Original full-master figures, marker/source/diversity and temporal analyses; source audit |
| `results/part_a/selected_original/` | Frozen 3,434-genome follow-up figures and statistics |
| `results/part_a/cohort/` | `Finalupdated_SEN_Genomes1_QC_merged.csv` with 3,434 rows, 69 original + 9 QC/provenance columns |
| `results/part_a/reconciled_metadata/` | Enriched full-master and selected+QC metadata; cross-population reconciliation |
| `results/part_a/temporal_top5/` | Master vs selected, seven mutually exclusive genotype categories by temporal block |
| `results/part_a/full_master_original/international_animal_source/` | International N=71 decomposed-source counts/percentages and detailed Excel audit |
| `results/part_a/master_33831/original_execution/` | Extracted original master/selected scripts and run/validation logs |

## Reproducibility, audit and limitations

- `scripts/check_original_rmd_contract.py` enforces the uploaded original's
  SHA-256 and full-master/selected chunk boundaries.
- `scripts/validate_full_master_repro.py` independently validates the master
  results, frozen selected cohort and QC handoff.
- `scripts/validate_enriched_metadata_outputs.py` validates matching keys
  and additive metadata/temporal-output consistency.
- The reproduced full-master animal-source breakdown may differ from a
  previously reported manuscript Table S1 (notably the Other-US Swine/Others
  split). Such differences are **flagged [REVIEW REQUIRED]**, not silently
  rewritten or automatically adopted as revised biological conclusions.
- An archival notebook is preserved for provenance, and the generated analysis
  code is the run target. A successful static GitHub check does **not** imply
  that the 33,831-genome biological analysis has been rerun on AjuloPC.

**Scope boundary:** no `analysis/clades/verify_table_s5_*.py`,
`scripts/run_amrfinder_table_s5_audit.sh`, Table S5 count data, or Table S5
validation results are included in this Part A2-only update. This distinction
does not remove or overwrite that separate analysis on other branches.
