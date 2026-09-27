# Publication audit and freeze gate

This document is the publication-facing reconciliation record for
`SEN_genomics_pipeline`. The authoritative historical source snapshot is
`ajulojays/SENBio-Final@78abf1e0efb3ed2eb91d7de72fc928a3b90e64f5`.

The cleaned repository has passed an end-to-end four-isolate smoke test through
reads/QC, typing, core-SNP phylogeny, assembly/pangenome, and
AMR/virulence/plasmid screening. Smoke validation demonstrates runtime
reproducibility; it does **not** prove that the full-study results were
historically generated with every currently pinned package version.

## 1. Cohort and metadata

**Status: RESOLVED**

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

`workflow/05_phylogeny/20_tree_clades.sh` is now a validation stage. It no
longer re-derives clades from historical `Clade1.txt/Clade2.txt/Clade3.txt`
lists. It requires the IQ-TREE SRR set to match the canonical final clade
metadata exactly.

## 2. Read/QC and typing

**Status: RESOLVED for reproducible reconstruction**

The smoke pipeline validates:

- SRA retrieval
- fastp trimming
- FastQC/MultiQC
- Kraken2 filtering to Enterobacteriaceae taxid 543 + descendants
- post-filter depth >=30x
- SeqSero2
- SKESA/SISTR
- MLST

SISTR compatibility is pinned to Python 3.10 with a setuptools release that
still provides `pkg_resources`.

### Historical version question

**OPEN:** manuscript MLST version = 2.23.0; reconstructed environment =
2.33.1. No source-code evidence in the frozen repository establishes which
version produced the historical final MLST outputs. Do not silently replace one
with the other in the manuscript.

## 3. Core SNP calling

**Status: RESOLVED for command contract**

Reference: P125109 / NC_011294.1.

Snippy 4.6.0 is now called with manuscript thresholds explicitly rather than
relying on defaults:

- `--mincov 10`
- `--minqual 100`
- `--mapqual 60`
- `--basequal 13`
- `--minfrac 0` (Snippy AUTO behavior)

The reconstructed alignment toolchain is pinned to BWA 0.7.18, SAMtools 1.20
and FreeBayes 1.3.6.

SnpEff remains a separate downstream annotation stage and must be validated
against the manuscript/reference database contract before final freeze.

## 4. Recombination filtering and final phylogeny

**Status: PARTIALLY RESOLVED**

The smoke pipeline validates Gubbins -> final exclusion step -> IQ-TREE 2.2.6,
1000 UFBoot2, BNNI, and ASC-aware model search.

The four historical post-Gubbins exclusions remain frozen in
`config/final_tree_exclusions.txt`.

### Historical version question

**OPEN:** manuscript Gubbins version = 3.0; reconstructed environment = 3.4.3.
The frozen source repository does not currently provide sufficient provenance to
decide which version produced the historical final recombination-filtered
alignment.

The manuscript reports TVM+F+ASC+R2 as the selected final IQ-TREE model. The
workflow correctly invokes ModelFinder with `MFP+ASC`; the selected model must
be verified from the historical/full-study IQ-TREE report before publication
freeze.

## 5. Assembly and pangenome

**Status: RESOLVED for workflow contract**

Publication path:

Shovill/SPAdes -> Prokka (*Salmonella enterica*) -> Panaroo strict -> MAFFT.

Definitions encoded in the analysis:

- accessory genes: 1-99% frequency
- statistical loci must occur in at least 10 genomes
- universally present loci are not tested
- clade-defining accessory genes: >=90% within clade and <=5% outside
- FDR: Benjamini-Hochberg

Runtime package installation has been removed from the publication-facing
Panaroo-by-clade script.

## 6. AMR, virulence and plasmids

**Status: RESOLVED for reproducible reconstruction**

Presence threshold:

- coverage >=80%
- identity >=90%

AMRFinderPlus:

- software: 4.2.5
- database: 2026-08-07.1
- `--plus --organism Salmonella`
- exact database release is frozen; analysis never calls a runtime updater

ABRicate 1.2.0 database inventory:

- ResFinder: 3,206 sequences, 2025-Dec-5
- VFDB: 4,592 sequences, 2025-Dec-5
- PlasmidFinder: 488 sequences, 2025-Dec-5

### Historical version question

**OPEN:** manuscript states AMRFinderPlus v3, whereas the reproducible
reconstruction uses 4.2.5. The smoke test establishes the current reconstruction
only; it does not establish that the historical full-study output was produced
with 4.2.5.

## 7. Clade-stratified statistics

**Status: CONTRACT CORRECTED; FULL-STUDY RERUN REQUIRED**

Publication-facing AMR, ResFinder, VFDB and plasmid scripts now enforce:

- exclude loci present in <10 genomes
- do not test zero/full non-variable strata
- omnibus chi-square
- Monte Carlo p-value when any expected cell is <5
- fixed 5,000 Monte Carlo resamples with RNG seed 12345
- Benjamini-Hochberg FDR across tested loci
- pairwise two-sided Fisher exact tests after significant omnibus result
- Holm adjustment within locus

The earlier code attempted an unsupported `chi2_contingency(...,
method="montecarlo", num_resamples=5000)` call and silently fell back to the
asymptotic chi-square result. That has been corrected to SciPy's
`MonteCarloMethod` API. Because this can alter sparse-count p-values, the
full-study statistical tables must be regenerated before publication freeze.

## 8. SNP clade enrichment

**Status: CONTRACT CORRECTED; FULL-STUDY RERUN REQUIRED**

Encoded manuscript rules:

- FDR <=0.05
- clade-defining threshold: >=40% within clade and <1% outside
- bins: 40-<60%, 60-<80%, 80-<90%, >=90%
- loci present in <10 genomes are excluded from inferential testing
- loci present in all genomes are excluded from inferential testing
- canonical split-clade metadata is the default input

Both combined-Clade-1 and split-Clade-1A/1B analyses may be retained as
sensitivity/descriptive outputs, but manuscript claims using the final clade
scheme must use the split metadata.

## 9. Remaining publication freeze gates

The branch should **not** be merged to `main` as a final publication freeze
until all of the following are complete:

1. rerun full-study clade AMR/ResFinder/VFDB/plasmid statistics after the
   Monte Carlo correction;
2. rerun full-study SNP clade statistics with the >=10/non-universal locus
   filter;
3. rerun or validate full-study Panaroo-by-clade outputs under the same locus
   filter;
4. validate SnpEff annotation/database provenance;
5. verify the full-study IQ-TREE report selected TVM+F+ASC+R2;
6. resolve, or explicitly document in the manuscript, historical version
   discrepancies for MLST, Gubbins and AMRFinderPlus;
7. run `bash scripts/publication_audit.sh` and require PASS;
8. review full-study regenerated tables against manuscript numbers before merge.

## 10. Merge policy

A smoke PASS is necessary but not sufficient for publication freeze. The final
merge should occur only after the regenerated full-study statistical outputs are
checked against the manuscript and all remaining OPEN provenance items have
either been resolved with evidence or explicitly described as reconstructed
software provenance.
