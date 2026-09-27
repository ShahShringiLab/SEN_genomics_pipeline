# Fresh-Linux smoke test

This guide reproduces the **four-isolate runtime smoke test** from a new Linux
machine. The smoke cohort is intentionally small and is used to validate
software installation, database bootstrapping, file handoffs, checkpoint/resume
behavior, and downstream analysis execution. It is **not** used for biological
inference.

## 1. Hardware and disk

The four genomes themselves are small, but the smoke test intentionally uses
the same pinned Kraken2 Standard database as the manuscript workflow.

Plan for:

- Linux x86_64 (Ubuntu 22.04/24.04 or comparable)
- >=16 GB RAM; 32 GB or more is preferable
- >=8 CPU threads; more will reduce runtime
- >=200 GB free disk for first-time Kraken2 download/extraction plus working
  files

The pinned Kraken2 database is approximately 80 GB compressed and approximately
103 GB extracted. The installer removes the archive after successful validation
unless configured otherwise.

## 2. Install system prerequisites

On Ubuntu/Debian:

```bash
sudo apt update
sudo apt install -y git curl wget ca-certificates bzip2 build-essential
```

## 3. Install Miniforge

The repository runners use `conda`.

```bash
cd ~
curl -L -o Miniforge3.sh \
  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

bash Miniforge3.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init bash
```

Open a new shell, or run:

```bash
source "$HOME/miniforge3/etc/profile.d/conda.sh"
```

Confirm:

```bash
conda --version
git --version
```

## 4. Clone the manuscript repository

```bash
git clone https://github.com/ShahShringiLab/SEN_genomics_pipeline.git
cd SEN_genomics_pipeline
```

If reproducing a specific manuscript release, checkout the tag/commit stated in
the manuscript instead of using a moving branch.

## 5. Validate the repository

```bash
bash scripts/validate_repo.sh
python scripts/validate_analysis_contracts.py
python tests/smoke_test/validate_smoke.py
```

Expected cohort:

| sample | clade |
|---|---|
| SRR1033488 | Clade 1A |
| SRR1220773 | Clade 1B |
| SRR10007520 | Clade 2 |
| SRR10005236 | Clade 3 |

The P125109 reference, NC_011294.1, is bootstrapped separately.

## 6. Run the complete smoke test

The smoke runner creates required Conda environments as needed, downloads the
four SRA runs, bootstraps pinned databases/reference assets, and resumes from
validated checkpoints.

```bash
bash tests/smoke_test/run_smoke_pipeline.sh
```

The complete smoke path is:

```text
SRA download
-> raw FastQC/MultiQC
-> fastp
-> trimmed FastQC/MultiQC
-> coverage
-> Kraken2 + Enterobacteriaceae extraction
-> SeqSero2
-> SKESA + SISTR
-> MLST
-> P125109 reference bootstrap
-> Snippy + snippy-core
-> Gubbins
-> final alignment cleanup
-> IQ-TREE
-> Shovill
-> Prokka
-> Panaroo
-> AMRFinderPlus
-> ABRicate ResFinder/VFDB/PlasmidFinder
-> analysis/clades
-> Panaroo-by-clade
-> analysis/snps
```

The R files `analysis/snps/FGA.R` and `FGA2.R` are syntax-parsed during the
smoke test. Their full COG analyses require external EggNOG outputs and are
therefore part of the full manuscript workflow rather than four-genome
biological inference.

## 7. Expected final message

A successful run ends with:

```text
SEN END-TO-END SMOKE PIPELINE: PASS
```

and the downstream analysis block reports:

```text
SMOKE DOWNSTREAM ANALYSIS RESULT: PASS
```

Because the smoke cohort has only four isolates, the manuscript rule excluding
loci present in fewer than 10 genomes means many inferential result tables will
be empty. That is expected.

## 8. Resume after interruption

Simply rerun:

```bash
bash tests/smoke_test/run_smoke_pipeline.sh
```

Validated completed stages are skipped using checkpoint signatures under:

```text
tests/smoke_test/work/.checkpoints/
```

To force every stage to rerun:

```bash
SEN_FORCE_RERUN=1 bash tests/smoke_test/run_smoke_pipeline.sh
```

To prevent adoption of valid pre-existing outputs:

```bash
SEN_ADOPT_EXISTING_OUTPUTS=0 bash tests/smoke_test/run_smoke_pipeline.sh
```

## 9. Run individual smoke blocks

```bash
bash tests/smoke_test/run_reads_qc.sh
bash tests/smoke_test/run_kraken_typing.sh
bash tests/smoke_test/run_core_phylogeny.sh
bash tests/smoke_test/run_assembly_pangenome.sh
bash tests/smoke_test/run_amr_vf_plasmid.sh
bash tests/smoke_test/run_analysis_layer.sh
```

## 10. Important interpretation

A smoke PASS means the workflow executes correctly on the representative
four-isolate cohort. It does **not** reproduce the manuscript's clade
frequencies, SNP counts, pangenome statistics, enrichment results, or final
3,306-isolate phylogeny.
