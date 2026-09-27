# Smoke milestone 1 — reads, QC, Kraken2 and typing

Status: **PASS**

Representative cohort: four *Salmonella Enteritidis* isolates, one from each
final major clade.

Validated stages:

1. SRA retrieval and paired FASTQ generation
2. raw FastQC + MultiQC
3. fastp trimming
4. trimmed FastQC + MultiQC
5. coverage estimation
6. Kraken2 classification
7. KrakenTools Enterobacteriaceae (taxid 543 + descendants) extraction
8. post-Kraken depth validation
9. SeqSero2
10. SKESA + SISTR
11. MLST

Observed smoke results:

| Sample | Initial depth | Post-Kraken depth | Purity |
|---|---:|---:|---:|
| SRR10007520 | 64.11x | 63.00x | 98.07% |
| SRR10005236 | 72.48x | 71.74x | 98.95% |
| SRR1033488 | 79.50x | 78.39x | 98.50% |
| SRR1220773 | 117.80x | 116.90x | 99.09% |

All four isolates:

- retained >=30x post-Kraken depth;
- were called Enteritidis by SeqSero2;
- were called Enteritidis by SISTR;
- produced SKESA assemblies;
- produced MLST results.

This milestone validates runtime contracts and reproducibility behavior only. It
does not constitute biological validation of the full 3,306-isolate study.

The smoke wrappers are checkpointed: validated completed stages are skipped,
while incomplete stages or stages whose script/environment signature changed
are rerun.
