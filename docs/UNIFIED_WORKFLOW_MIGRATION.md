# Unified Part A / Part B directory migration

**Status:** preparation branch only — NOT merged into `main`.

This branch is dedicated to consolidating analysis and workflow scripts into a single `workflow/` hierarchy, with **original filenames unchanged**.

| Current location | Intended location |
| --- | --- |
| `analysis/part_a/01_phenotype_trends.R` | `workflow/A1_Phenotype/01_phenotype_trends.R` |
| Other files in `analysis/part_a/` | `workflow/A2_Genotype/` (except the A1 phenotype script) |
| `workflow/01_reads/`, `workflow/02_qc/` and trimming under `workflow/03_typing/` | `workflow/B1_Reads_QC/` |
| Remaining `workflow/03_typing/` scripts | `workflow/B2_Typing/` |
| `workflow/04_core_snp/`, `workflow/05_phylogeny/`, `analysis/snps/` and clade-phylogeny files | `workflow/B3_Phylogenomics/` |
| `workflow/06_assembly_pangenome/` | `workflow/B4_Pangenome/` |
| `workflow/07_amr_vf_plasmid/` and clade AMR/plasmid/VFDB files | `workflow/B5_AMR_Virulence_Plasmids/` |

**Protected paths:** `data/`, `metadata/`, `provenance/`, `reference/`, `environments/` and all original filename basenames.

### Required migration checks

1. Inventory every executable and its path references, including the `scripts/` launchers, smoke tests, configuration, README and supplementary documentation.
2. Use `git mv` to retain version history and preserve original filenames. Move `analysis/part_a/lib/common.R` into `workflow/A2_Genotype/lib/common.R`; retain or explicitly update `workflow/lib/common.sh` references.
3. Rewrite **path references**, not scientific logic or output artifacts. Account for imports and scripts that infer `ROOT` relative to `__file__` or their working directory.
4. Confirm no old executable paths remain in launchers/tests. Check all referenced new files exist and all shell scripts remain executable.
5. Run `bash scripts/publication_audit.sh`, `bash scripts/validate_repo.sh`, the canonical Part A runtime checks, and the full Part B smoke test in a compatible test environment.
6. Compare result counts and designated outputs to the pre-migration reference. Only merge after all checks pass; if full runtime testing is unavailable, leave the changes on this branch.

The existing README is the authoritative documentation structure: retain all Part A1, Part A2, and Part B1–B5 descriptions while changing only paths and repository-layout text.
