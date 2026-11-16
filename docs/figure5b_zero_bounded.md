# Recreate Figure 5B without negative violin tails

This branch packages the final SEN IQ-TREE Newick tree supplied for the Figure 5B correction. It retains the canonical clade assignments in the repository and uses the original `r_snps.R` approach.

## Quick start (Ubuntu / WSL / macOS with Conda)

```bash
git clone --branch figure5b-zero-bounded-violin --single-branch https://github.com/ShahShringiLab/SEN_genomics_pipeline.git
cd SEN_genomics_pipeline
bash scripts/run_figure5b.sh
```

The runner creates `sen_figure5b` from `environments/figure5b.yaml` if it is not already installed. It then runs the R script.

## Inputs (included in the branch)

- `data/phylogeny/SSLAB_FINAL.treefile`: the supplied 3,307-tip IQ-TREE Newick tree (3,306 SRR records and the Reference)
- `metadata/final_clade_metadata.tsv`: canonical clade assignments (including Clade 1A, Clade 1B, Clade 2, Clade 3)
- `metadata/SEN_Genomes.csv`: used for the original master-metadata sample restriction

No manual input-file copying is required.

## Outputs (ignored by Git)

`results/figure5b_zero_bounded/`:

- `Figure5B_violin_zero_bounded.tif` (300 DPI, LZW)
- `Figure5B_violin_zero_bounded.pdf` (vector)
- `Figure5B_violin_zero_bounded.png` (300 DPI preview)
- `Figure5B_distance_summary.csv` (genomes, pairs, min, median, max)

## Specific correction

The historical script used `geom_violin(trim = FALSE)`, allowing the Gaussian KDE to extend beyond observed values. Here the exact violin-density input is preserved, but `trim = TRUE` and a nonnegative display boundary prevent tails below zero. The lineage order and colors match the original Figure 5B.

Distances are **tree-derived patristic branch-length sums × 67,736** (the recombination-masked SNP-only alignment length), exactly as the historical `r_snps.R` plot calculated them. They are not independently counted pairwise SNP differences.

Historical source: [original `r_snps.R`](https://github.com/ajulojays/SENBio-Final/blob/78abf1e0efb3ed2eb91d7de72fc928a3b90e64f5/r_snps.R#L222-L355).

## Custom paths

For another tree, clade table, or output path:

```bash
bash scripts/run_figure5b.sh \
  /path/to/SSLAB_FINAL.treefile \
  metadata/final_clade_metadata.tsv \
  results/figure5b_zero_bounded \
  metadata/SEN_Genomes.csv
```
