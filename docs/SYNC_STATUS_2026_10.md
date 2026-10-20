# October 2026 SEN manuscript synchronization — staging review

This branch is a **draft publication candidate**, not a validated release. It retains the lab's original Part B pipeline and historical software-provenance audit while integrating the private manuscript-release Part A implementation and the latest clade-assignment developments.

## Exact provenance

- Public lab base: `ShahShringiLab/SEN_genomics_pipeline@63989a22b09f8b8958bc47d36350b08df0c7c184` (2026-09-27).
- Part A manuscript source: `ajulojays/SEN_genomics_pipeline@d2f34ca4e28918a793c302a2dcab9ecce62949b7` (`manuscript-release`, 2026-10-03).
- Latest clade-script/sidecar source: `ajulojays/SEN_genomics_pipeline@2cc69509f35cb34c1ce51d366f83fac3c2365394` (`main`, 2026-10-03).

## Cohort and scientific contracts

- Master cohort 33,831; selected cohort 3,434 (Run+Isolate keys).
- Final phylogeny: 3,306 study isolates + P125109 reference, with 128 selected isolates outside the final study phylogeny.
- Corrected QR/FQR core marker status: 2,288 marker-positive, 1,146 marker-negative.
- Three historical QR status discrepancies are **preserved in an explicit audit**, not silently rewritten in the archived source cohort.
- Curated QR/clade table uses a **blank** clade for the 128 excluded isolates. The complementary compact clade sidecar labels them **Not in phylogeny**. Do not treat this as a sixth evolutionary clade.

## Deliberate reconciliation

- Adopt `manuscript-release` Part A logic, based on curated `AMR genotypes core`, **not** the development `main` branch's broader `AMR genotypes` token calls.
- Retain the validated October 3 improvements to `analysis/clades/14b_cladesummary.py`.
- Port the development fix for `dplyr::if_else` group-size recycling in `05_diversity_composition.R`.
- Port the `NA`-to-empty selection-key normalization fix into the manuscript input validator.
- Keep the lab's more detailed `docs/version_audit.md` for historical software-version provenance. Observed workstation versions are not automatically identical to historical execution versions.
- Preserve canonical reference, clade metadata, existing Part B workflow, and four-genome smoke-test code.

## Blocking release inputs — do not merge yet

Two bulk CSVs have **not** been transferred to this public staging branch because the connected API bulk transfer failed. They remain in the private manuscript-release repository:

- `data/part_a/genotype/US-SEN_12192025.csv` — 33,831 records; Git blob `6ffc47814dca722ac75e7b335d4b976dbfaad9ee` (34,087,900 bytes).
- `data/part_a/selection/Finalupdated_SEN_Genomes1.csv` — 3,434 records; Git blob `b51383acf25d52fbd2eb979a416092acd8d2d2ca` (4,197,893 bytes).

**Release requirement:** Confirm that these derived source records are authorized for public publication, transfer intact through a Git client, verify matching blob hashes, then run the relevant checks. The committed smaller phenotype, selection-key, and QR/clade audit inputs are not substitutes for the missing files.

## Validation required before public main

```bash
bash scripts/validate_repo.sh
python3 scripts/validate_part_a_inputs.py
bash scripts/run_part_a.sh
bash scripts/publication_audit.sh
bash tests/smoke_test/run_smoke_pipeline.sh
```

The smoke workflow requires its declared databases/software and may involve read retrieval. The full manuscript run is documented separately in `docs/FULL_RUN_FRESH_LINUX.md`.

No `SEN PART A: PASS`, `STATIC PUBLICATION AUDIT: PASS`, or end-to-end smoke-test `PASS` has been asserted for this staging synchronization. Its status stays **draft/pending** until required inputs and the validation evidence are available.
