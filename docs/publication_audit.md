# Publication release audit

This document records the publication-facing validation state of
`SEN_genomics_pipeline`. The repository contains the complete SEN workflow
implementation, with pinned analysis environments, canonical metadata,
database/reference bootstrap utilities, smoke tests, and full-run instructions.

## 1. Part A phenotype/genotype/cohort workflow

**Status: WORKFLOW IMPLEMENTED**

Part A production scripts are under `analysis/part_a/` and are executed with:

```bash
bash scripts/run_part_a.sh
```

A completed release requires the phenotype input, upstream genotype/metadata table, final 3,434-isolate selection key, and committed final 3,434-genome cohort under `data/part_a/`. The bulk genotype and final-cohort CSVs are not yet present on this staging branch.

Runtime Part A validation requires:

- phenotype shared-denominator consistency before trend calculation;
- exactly 3,434 unique isolates in the selection key;
- reconstruction of the selected cohort from the upstream genotype table;
- exact agreement of rebuilt collection year and source/origin classification with the committed final cohort;
- explicit audit of the three known notebook-rule versus canonical QR_genotype differences;
- use of curated AMRFinderPlus-core-derived corrected QR status (2,288 positive, 1,146 negative) for downstream production summaries, preserving the three historical discrepancies;
- deterministic bootstrap/permutation seeds for diversity analyses;
- post-run validation of cohort size/uniqueness, mismatch counts, source totals,
  phenotype result dimensions, diversity source coverage, and required figures
  before the runner reports PASS.

The original exploratory notebook is preserved under `provenance/SEN_Genomes2.Rmd`; the modular R scripts are the executable release implementation. The archived master genotype CSV is also preserved unchanged; known quoting/escaped-comma artifacts are normalized only in memory for parsing.

## 2. Cohort and metadata

**Status: VERIFIED**

Canonical study metadata:

- `metadata/SEN_Genomes.csv`: 3,434 selected genomes
- `metadata/final_clade_metadata.tsv`: 3,306 study SRRs + 1 reference
- final clade counts:
  - Clade 1A: 169
  - Clade 1B: 759
  - Clade 2: 894
  - Clade 3: 1,357
  - Unassigned: 127
  - Reference: 1

`workflow/05_phylogeny/20_tree_clades.sh` validates that the IQ-TREE tip set
matches the canonical final clade metadata.

## 3. Read/QC and typing

**Status: WORKFLOW VERIFIED**

The release workflow includes:

- SRA retrieval
- fastp trimming
- FastQC/MultiQC
- Kraken2 filtering to Enterobacteriaceae taxid 543 + descendants
- post-filter depth >=30x
- SeqSero2
- SKESA/SISTR
- MLST

Software versions for the release are defined in `environments/` and summarized
in `docs/version_audit.md`.

## 4. Core SNP calling and annotation

**Status: WORKFLOW VERIFIED**

Reference: P125109 / NC_011294.1.

The release uses Snippy 4.6.0 with explicit thresholds:

- `--mincov 10`
- `--minqual 100`
- `--mapqual 60`
- `--basequal 13`
- `--minfrac 0` (Snippy AUTO)

The release environment also pins BWA 0.7.18, SAMtools 1.20, FreeBayes 1.3.6,
and SnpEff 5.1. The SnpEff database is built from the bootstrapped P125109
GenBank record.

## 5. Recombination filtering and phylogeny

**Status: WORKFLOW VERIFIED**

The release implements:

- Gubbins recombination filtering
- final alignment exclusions from `config/final_tree_exclusions.txt`
- IQ-TREE 2.2.6
- ascertainment-bias-aware model selection
- 1,000 UFBoot2 replicates
- BNNI optimization

The manuscript reports `TVM+F+ASC+R2` as the selected full-study model; the
full-run guide includes verification of this model from the generated IQ-TREE
report.

## 6. Assembly and pangenome

**Status: WORKFLOW VERIFIED**

Publication path:

Shovill/SPAdes -> Prokka (*Salmonella enterica*) -> Panaroo strict -> MAFFT.

Analysis definitions:

- accessory genes: 1-99% frequency
- inferential loci must occur in at least 10 genomes
- universally present loci are excluded from inferential testing
- clade-defining accessory genes: >=90% within clade and <=5% outside
- Benjamini-Hochberg FDR correction

## 7. AMR, virulence and plasmids

**Status: WORKFLOW VERIFIED**

Presence thresholds:

- coverage >=80%
- identity >=90%

The release uses:

- NCBI AMRFinderPlus 4.2.5
- AMRFinderPlus database 2026-08-07.1
- `--plus --organism Salmonella`
- ABRicate 1.2.0
- ResFinder, VFDB and PlasmidFinder database snapshots recorded by the database
  bootstrap workflow

Database updates are separated from analysis execution so a run uses the
configured release assets.

## 8. Clade-stratified statistics

**Status: ANALYSIS CONTRACT VERIFIED**

Publication-facing scripts enforce:

- exclude loci present in <10 genomes
- exclude zero/full non-variable loci from inferential testing
- omnibus chi-square
- Monte Carlo p-value when expected cells are <5
- 5,000 Monte Carlo resamples with RNG seed 12345
- Benjamini-Hochberg FDR across tested loci
- pairwise two-sided Fisher exact tests after a significant omnibus result
- Holm adjustment within locus

## 9. SNP clade enrichment

**Status: ANALYSIS CONTRACT VERIFIED**

Encoded rules:

- FDR <=0.05
- within-clade presence >=40%
- outside-clade presence <1%
- bins: 40-<60%, 60-<80%, 80-<90%, >=90%
- loci present in <10 genomes excluded from inferential testing
- universal loci excluded from inferential testing
- canonical split-clade metadata used by default

## 10. Release verification

For manuscript-scale verification:

1. run the full workflow using `docs/FULL_RUN_FRESH_LINUX.md`;
2. compare regenerated full-study tables with manuscript values;
3. verify the selected IQ-TREE model in the generated report;
4. retain SnpEff and external EggNOG/PHASTER/HMMER provenance with the outputs;
5. run `bash scripts/publication_audit.sh` and require PASS.

## 11. Release model

The ShahShringiLab repository represents the complete manuscript-facing SEN
workflow. Development and testing are maintained separately; validated updates
are promoted into the manuscript repository as complete release snapshots.
