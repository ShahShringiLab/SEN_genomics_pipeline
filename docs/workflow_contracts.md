# Workflow contracts

This document defines the publication-facing handoff between stages. Paths are defaults; local/HPC paths can be overridden through `config/paths.env` or the corresponding `SEN_*` environment variables.

## Part A1. Phenotype resistance trends

**Script:** `analysis/part_a/01_phenotype_trends.R`  
**Environment:** `environments/12_stats.yaml`  
**Input:** `data/part_a/phenotype/FDASEN_AMRPhenotype2.csv`  
**Output:** phenotype trend table and Human-versus-Chicken figures under `results/part_a/phenotype/`  
**Key contract:** 2024 and 2025 are excluded from the plotted analysis; the historical shared denominator is used only after verifying that Ciprofloxacin-tested and DSC-tested counts are identical for each record; normal-approximation 95% confidence intervals are reproduced.

## Part A2. Genotype epidemiology and final cohort

**Scripts:** `analysis/part_a/02_build_final_cohort.R`, `03_genotype_marker_summary.R`, `04_source_profiles.R`, `05_diversity_composition.R`  
**Environment:** `environments/12_stats.yaml`  
**Inputs:** master genotype/metadata table, final 3,434-isolate selection key, and committed final 3,434-genome cohort under `data/part_a/`  
**Outputs:** rebuilt/validated final cohort, QR/FQR marker summaries, temporal/source tables, source/origin profiles, Hill-number diversity estimates, bootstrap confidence intervals and permutation tests under `results/part_a/`  
**Key contract:** the final selection key contains 3,434 unique isolates. Collection year and source/origin classification rebuilt from upstream data must agree exactly with the committed final cohort. The notebook-derived broad QR/FQR marker rule has three documented differences from the committed final `QR_genotype`; these are written to an audit table and the committed final QR status is used for binary downstream summaries, while marker-level summaries remain token-based. The archived master CSV remains unchanged and parser normalization is performed only in memory. `scripts/validate_part_a_outputs.py` must pass before the Part A runner reports PASS. The final Part A cohort is the starting cohort for Part B.

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
by the reproducible Kraken environment.  
**Key contract:** TaxID 543 is retained; downstream study inclusion requires
post-filter depth >=30x.

## 4. Species confirmation and typing

**Scripts:** `10_seqsero2.sh`, `11_sistr.sh`, `12_mlst.sh`  
**Environments:** `03_seqsero2.yaml`, `03_sistr.yaml`, `04_mlst.yaml`  
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
are depth 10, VCF QUAL 100, MAPQ 60 and base quality 13. Snippy AUTO allele
fraction behavior is invoked with `--minfrac 0`.

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
**Key contract:** AMRFinderPlus runs with `--plus --organism Salmonella`
against an explicit frozen database directory; the analysis stage never performs
a database update. Database bootstrap is separated into
`scripts/setup_amr_databases.sh`, which captures the resolved AMRFinderPlus
database release and full ABRicate database inventory as provenance. ABRicate
runs ResFinder, VFDB and PlasmidFinder with the manuscript binary-presence
thresholds encoded directly at screening time: coverage >=80% and identity
>=90%. ABRicate database updates are not performed at runtime; the databases
bundled with the pinned ABRicate package are inventoried and validated before
screening.

## 10. Downstream clade/SNP statistics

**Scripts:** `analysis/clades/*`, `analysis/snps/*`  
**Environment:** `12_stats.yaml`  
**Input:** canonical clade assignments and stage-specific SNP/gene/element
tables  
**Output:** clade summaries, enrichment tables, source comparisons and figures  
**Key manuscript contract:** chi-square omnibus testing, Monte Carlo treatment
when expected counts are low, pairwise Fisher exact tests with Holm correction,
FDR across loci, and exclusion of loci seen in <10 genomes or universally
present loci.



### SISTR CLI migration correction

The publication-facing SISTR stage uses `-i <fasta> <genome_name>` for explicit
genome naming. Historical code used `-n <sample>`, but in `sistr_cmd`
`-n/--novel-alleles` specifies a novel-alleles FASTA output path rather than a
genome name. The corrected stage also preserves completed SKESA assemblies,
reruns only missing SISTR results, and retains per-sample SISTR logs on failure.
