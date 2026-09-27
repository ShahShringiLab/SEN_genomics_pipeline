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


## Download performance

Large HTTP/S3 assets use `aria2c` when available, with resumable segmented
downloads (16 connections by default), and fall back to `wget -c` only when
necessary. Override the segment count with `SEN_ARIA2_CONNECTIONS`.

SRA acquisition uses parallel `prefetch` plus parallel `fasterq-dump`, with
concurrency derived conservatively from host CPU count. Users can still override
`SEN_PREFETCH_JOBS`, `SEN_DUMP_JOBS`, and `SEN_THREADS_PER_DUMP` for a
specific workstation or HPC node.


## Intelligent resume and checkpoints

Smoke stages are idempotent and resume-aware. Before running an expensive
stage, the wrapper validates the complete expected output contract for every
smoke sample. A stage is skipped only when those outputs are valid and its
stored workflow/config signature matches the current scripts and environment
definition.

Validated outputs produced before checkpointing was introduced are adopted once
by default, so an interrupted pilot can continue without recomputing completed
work. If a scientific script or environment YAML changes later, the signature
changes and that stage reruns automatically.

Controls:

```bash
# Force every stage to rerun even when outputs/checkpoints are valid.
SEN_FORCE_RERUN=1 bash tests/smoke_test/run_smoke_pipeline.sh

# Require explicit checkpoints; do not adopt pre-existing outputs.
SEN_ADOPT_EXISTING_OUTPUTS=0 bash tests/smoke_test/run_smoke_pipeline.sh
```

Checkpoint records live under
`tests/smoke_test/work/.checkpoints/` by default.


## Core-SNP and phylogeny block

The next checkpointed smoke block is:

```text
P125109 / NC_011294.1 reference bootstrap
→ Snippy per isolate
→ snippy-core
→ Gubbins recombination filtering
→ historical final-tree exclusion step
→ IQ-TREE 2.2.6 with ModelFinder+ASC, 1000 UFBoot2, BNNI
```

Run it alone with:

```bash
bash tests/smoke_test/run_core_phylogeny.sh
```

or run the currently implemented workflow from the beginning with:

```bash
bash tests/smoke_test/run_smoke_pipeline.sh
```

The reference bootstrap downloads the exact chromosome accession
`NC_011294.1` from NCBI and validates the downloaded sequence before use.
SnpEff functional annotation remains a separate downstream stage rather than
being coupled to Snippy/core-alignment construction.


## End-to-end smoke status

All currently publication-facing computational blocks now pass on the
four-isolate smoke cohort:

```text
Reads/QC
→ Kraken2
→ SeqSero2 + SISTR + MLST
→ Snippy + snippy-core
→ Gubbins
→ IQ-TREE
→ Shovill
→ Prokka
→ Panaroo
→ AMRFinderPlus
→ ABRicate ResFinder/VFDB/PlasmidFinder
```

Frozen AMR database contracts:

- AMRFinderPlus database: `2026-08-07.1`
- ABRicate ResFinder: 3,206 sequences, 2025-Dec-5
- ABRicate VFDB: 4,592 sequences, 2025-Dec-5
- ABRicate PlasmidFinder: 488 sequences, 2025-Dec-5

The one-command smoke entry point remains:

```bash
bash tests/smoke_test/run_smoke_pipeline.sh
```
