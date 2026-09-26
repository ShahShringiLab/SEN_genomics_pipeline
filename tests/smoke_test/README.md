# SEN pipeline smoke test

This smoke cohort is a **runtime validation set**, not a biological subsample for
inference. It contains one confirmed isolate from each major final clade plus
the P125109 reference.

## Cohort

| Sample | Final clade | Source group | QR group |
|---|---|---|---|
| SRR1033488 | Clade 1A | Other-US | QR |
| SRR1220773 | Clade 1B | Human | QR |
| SRR10007520 | Clade 2 | Human | NON-QR |
| SRR10005236 | Clade 3 | Chicken-US | QR |
| Reference | Reference | — | — |

All four isolate records in the canonical study metadata are Illumina,
paired-end runs and have existing assemblies. The cohort was chosen to exercise
all four major clade labels while keeping the test small.

## Purpose

The smoke test should verify stage contracts and runtime behavior:

1. read retrieval
2. trimming/QC
3. Kraken2 filtering
4. SeqSero2 and SISTR
5. MLST
6. Snippy/core alignment
7. Gubbins
8. IQ-TREE
9. Shovill
10. Prokka
11. Panaroo
12. AMRFinderPlus and ABRicate

The smoke test does **not** validate full-study biological results, clade
frequencies, enrichment statistics, or the final 3,306-isolate phylogeny.

## Use

Validate the cohort metadata first:

```bash
python tests/smoke_test/validate_smoke.py
```

For read download, point the workflow at:

```bash
export SEN_SRR_LIST="$PWD/tests/smoke_test/smoke_srrs.txt"
```

Run computational stages in their documented environments. Heavy steps can be
executed on the lab workstation/HPC even if repository preparation is done on a
laptop.


## One-command staged runners

After validating the manifest, the current smoke-test sequence is:

```bash
bash tests/smoke_test/run_reads_qc.sh
```

Then configure a Kraken2 database and run:

```bash
export SEN_KRAKEN_DB=/path/to/kraken_db
bash tests/smoke_test/run_kraken_typing.sh
```

The second runner validates KrakenTools extraction, post-filter depth >=30x,
SeqSero2, SKESA/SISTR dual Enteritidis typing, and MLST for all four smoke
isolates. It uses project-local temporary storage and laptop-safe parallelism by
default.


## Fresh-machine behavior

The smoke workflow is being designed so that a fresh Linux/WSL machine can
clone the repository and run the pipeline without manually hunting for
databases. External scientific assets are pinned and bootstrapped by repository
scripts.

For Kraken2, if no complete `SEN_KRAKEN_DB` is already configured,
`run_kraken_typing.sh` automatically downloads and verifies the pinned
Kraken2 Standard snapshot defined in `config/database_sources.env`.

The full Standard database is intentionally used by default rather than a
reduced MiniKraken database because the goal is publication-facing
reproducibility. The pinned 2026-06-26 archive is approximately 80 GB compressed
and the extracted index is approximately 103 GB, so first-time setup requires
substantial disk space and download time. Subsequent runs reuse the installed
database.

The current one-command entry point is:

```bash
bash tests/smoke_test/run_smoke_pipeline.sh
```

As the remaining reference/SnpEff/AMR database contracts are finalized, their
bootstrap stages will be added to this same entry point rather than requiring
manual setup.
