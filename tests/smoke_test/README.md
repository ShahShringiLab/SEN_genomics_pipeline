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
