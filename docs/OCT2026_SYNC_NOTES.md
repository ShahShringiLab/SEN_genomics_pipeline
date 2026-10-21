# October 2026 SEN synchronized review branch

This document tracks a **staging candidate**, not a finished public release.

## Source versions

- Public base: ShahShringiLab/SEN_genomics_pipeline at 63989a22b09f8b8958bc47d36350b08df0c7c184
- Private scientific manuscript source: ajulojays/SEN_genomics_pipeline branch manuscript-release (d2f34ca4e28918a793c302a2dcab9ecce62949b7)
- Latest clade analysis and sidecar: ajulojays/SEN_genomics_pipeline main (2cc69509f35cb34c1ce51d366f83fac3c2365394)

## Outstanding two large files

The GitHub connector cannot transfer these blobs directly. They must be copied
from the private `manuscript-release` source via a credentialed local Git client,
before running Part A or merging the candidate:

| Path | Original Git blob SHA |
| --- | --- |
| `data/part_a/genotype/US-SEN_12192025.csv` | `6ffc47814dca722ac75e7b335d4b976dbfaad9ee` |
| `data/part_a/selection/Finalupdated_SEN_Genomes1.csv` | `b51383acf25d52fbd2eb979a416092acd8d2d2ca` |

From a working clone of ShahShringiLab, on this review branch:

```bash
git remote add development https://github.com/ajulojays/SEN_genomics_pipeline.git
git fetch development manuscript-release
git switch sync-sen-oct2026-review
git restore --source=development/manuscript-release -- \
  data/part_a/genotype/US-SEN_12192025.csv \
  data/part_a/selection/Finalupdated_SEN_Genomes1.csv

git hash-object data/part_a/genotype/US-SEN_12192025.csv
git hash-object data/part_a/selection/Finalupdated_SEN_Genomes1.csv
# Compare the reported SHA values to the table above before proceeding.
git add data/part_a/genotype/US-SEN_12192025.csv \
  data/part_a/selection/Finalupdated_SEN_Genomes1.csv
git commit -m "Restore verified canonical Part A data inputs"
git push origin sync-sen-oct2026-review
```

If the `development` remote is already present, do not add it again.

## Scientific reconciliation

The full corrected QR/clade table and the clade sidecar match all 3,434
Run/Isolate/clade assignments. The three gyrA_D87Y-positive Clade 3 genomes
(SRR15338589, SRR16764591, SRR16767069) are marker-positive under the
AMRFinderPlus-core-derived Part A definition but receive curated QS/FQS
overrides in the latest clade summary. These are not interchangeable labels.
Freeze a single operational definition for each reported phenotype, genotype,
and figure before publication.

## Release gates

```bash
bash scripts/publication_audit.sh
bash scripts/run_part_a.sh
bash tests/smoke_test/run_smoke_pipeline.sh
```

Retain the historical software provenance warnings in
`docs/version_audit.md` and `docs/publication_audit.md`. Do not report a
completed smoke test, regenerated full-study statistics or a final validated
release solely on the basis of passing static source-level checks.
