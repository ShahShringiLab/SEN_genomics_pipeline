# Workflow contracts

This document defines the publication-facing handoff between stages. Paths are
defaults; local/HPC paths can be overridden through `config/paths.env` or the
corresponding `SEN_*` environment variables.

## 1. Read retrieval

**Script:** `workflow/01_reads/01_download_fastq.sh`  
**Environment:** `environments/01_reads_qc.yaml`  
**Input:** accession list (`SEN_SRR_LIST`)  
**Output:** paired compressed FASTQs under `SEN_RAW_READS`  
**Key contract:** each accession must yield `*_1.fastq.gz` and
`*_2.fastq.gz`.

## 2. Raw/trimmed QC and trimming

**Scripts:** `02_fastqc_raw.sh`, `03_fastqc_trimmed.sh`,
`09_trim_fastp.sh`  
**Environment:** `environments/01_reads_qc.yaml`  
**Input:** raw paired FASTQs  
**Output:** trimmed paired FASTQs plus FastQC/MultiQC and fastp reports  
**Key trimming parameters:** sliding window 4; mean Q20; minimum length 50;
paired-end correction.

## 3. Kraken2 filtering and depth QC

**Scripts:** `04_estimate_coverage.sh`, `05_two_step_filter.sh`,
`06_kraken_cleanup.sh`, `07_presnippy_qc.py`, `08_qc_filter.py`  
**Environment:** `environments/02_kraken.yaml`  
**Input:** trimmed paired FASTQs and configured Kraken2 database  
**Output:** Enterobacteriaceae-retained paired reads under
`Kraken_cleanup/clean_trimmed_fastq` and QC tables  
**Extraction implementation:** KrakenTools `extract_kraken_reads.py` is supplied
by the reproducible Kraken environment; the missing historical local helper is
not required in the publication-facing pipeline.  
**Key contract:** TaxID 543 is retained; downstream study inclusion requires
post-filter depth >=30x. Historical QC scripts contain additional diagnostics;
the final manuscript inclusion rule must remain distinguishable from those
diagnostics.

## 4. Species confirmation and typing

**Scripts:** `10_seqsero2.sh`, `11_sistr.sh`, `12_mlst.sh`  
**Environments:** `03_serotyping.yaml`, `04_mlst.yaml`  
**Input:** trimmed reads for SeqSero2; Kraken-cleaned reads for SKESA/SISTR  
**Output:** SeqSero2 summary, SISTR assemblies/report, MLST report  
**Key contract:** study inclusion uses dual Enteritidis confirmation from
SeqSero2 and SKESA/SISTR before core-SNP analysis.

## 5. Core SNP calling and annotation

**Scripts:** `13_snippy.sh`, `14_snpeff.sh`, `15_snpeff_parse.py`,
`16_genome_audit.py`  
**Environment:** `05_snippy_snpeff.yaml`  
**Input:** confirmed isolate list, Kraken-cleaned paired reads, P125109
reference FASTA/GFF, SnpEff database  
**Output:** per-sample Snippy outputs, core alignment/VCF, annotated VCF and
variant report  
**Key contract:** P125109/NC_011294.1 is the reference. Publication thresholds
are depth 10, VCF QUAL 100, MAPQ 60 and base quality 13; major-allele behavior
must be audited against the actual historical Snippy invocation before final
Methods freeze.

## 6. Recombination filtering and phylogeny

**Scripts:** `17_gubbins.sh`, `18_gubbins_cleanup.sh`, `19_iqtree.sh`,
`20_tree_clades.sh`  
**Environments:** `06_gubbins.yaml`, `07_iqtree.yaml`  
**Input:** Snippy core alignment  
**Output:** recombination-filtered alignment, cleaned final alignment, ML tree,
bootstrap support and clade assignments  
**Key contract:** final IQ-TREE uses ASC-aware model selection, 1000 UFBoot2
replicates and BNNI. The manuscript reports TVM+F+ASC+R2 as the selected model.
The four historical post-Gubbins exclusions are tracked in
`config/final_tree_exclusions.txt`. Final clade metadata is
`metadata/final_clade_metadata.tsv`.

## 7. Assembly and annotation

**Scripts:** `21_shovill.sh`, `22_prokka.py`  
**Environments:** `08_shovill.yaml`, `09_prokka.yaml`  
**Input:** quality-controlled paired reads  
**Output:** Shovill draft assemblies followed by Prokka GFF annotations  
**Key contract:** the publication path is Shovill -> Prokka. Prokka is
configured for *Salmonella enterica*. Environment creation is external to the
analysis script.

## 8. Pangenome

**Scripts:** `23_panaroo.sh`, `23_panaroo_byclade.py`  
**Environment:** `10_panaroo.yaml`  
**Input:** Prokka GFF files  
**Output:** strict-mode Panaroo pangenome/core alignment and clade accessory-gene
analysis  
**Key contract:** MAFFT is used for core alignment. Publication definitions:
accessory genes are 1-99% frequency; clade-defining genes are >=90% within the
clade and <=5% outside.

## 9. AMR, virulence and plasmid screening

**Scripts:** `24_amrfinderplus.sh`, `25_abricate.sh`  
**Environment:** `11_amr_abricate.yaml`  
**Input:** draft assemblies  
**Output:** AMRFinderPlus reports and ABRicate ResFinder/VFDB/PlasmidFinder
reports  
**Key contract:** AMRFinderPlus runs with `--plus --organism Salmonella`.
Binary presence in manuscript analyses is coverage >=80% and identity >=90%.
The current AMRFinder workflow still performs `amrfinder -u`; this must be
replaced by a documented/pinned database snapshot before final publication
freeze.

## 10. Downstream clade/SNP statistics

**Scripts:** `analysis/clades/*`, `analysis/snps/*`  
**Environment:** `12_stats.yaml`  
**Input:** canonical clade assignments and stage-specific SNP/gene/element
tables  
**Output:** clade summaries, enrichment tables, source comparisons and figures  
**Key manuscript contract:** chi-square omnibus testing, Monte Carlo treatment
when expected counts are low, pairwise Fisher exact tests with Holm correction,
FDR across loci, and exclusion of loci seen in <10 genomes or universally
present loci. Each analysis script must be checked against this contract before
the final statistical freeze.

## Known unresolved provenance items

Current workstation package versions do not prove the exact historical
execution versions. The version discrepancies documented in
`docs/version_audit.md` remain open provenance questions, especially MLST,
Gubbins and AMRFinderPlus.

The smoke test validates execution contracts only. It is not expected to
reproduce full-study phylogenetic topology, clade statistics or enrichment
results from four isolates.


### SISTR CLI migration correction

The publication-facing SISTR stage uses `-i <fasta> <genome_name>` for explicit
genome naming. Historical code used `-n <sample>`, but in `sistr_cmd`
`-n/--novel-alleles` specifies a novel-alleles FASTA output path rather than a
genome name. The corrected stage also preserves completed SKESA assemblies,
reruns only missing SISTR results, and retains per-sample SISTR logs on failure.
