# Original Part A2: faithful Rmd preservation (33,831 genomes)

**Authoritative source:** `provenance/SEN_Genomes2_uploaded_original_20261021.Rmd` (the user's full 30-chunk source, SHA-locked and tracked in this branch). The earlier `provenance/SEN_Genomes2.Rmd` has 26 chunks and is **not interchangeable**. See [the current Part A2 release](PART_A2_UPDATED.md) for the canonical workflow. The complete original is included in this branch; no manual copy is needed.

```bash
git clone --branch feature/part-a2-updated-20261023 \
 https://github.com/ShahShringiLab/SEN_genomics_pipeline.git SEN_MasterA2
cd SEN_MasterA2
mkdir -p provenance
source ~/miniconda3/etc/profile.d/conda.sh
conda install -n sen_stats -c conda-forge r-writexl -y
bash scripts/run_master_part_a2.sh
```

The runner first normalizes CSV quoting/escaped commas and exports **exactly 33,831 rows × original 51 columns** as immutable input snapshot at
`results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.xlsx`
and equivalent CSV, with field-alignment audit and checksums.

Then `scripts/prepare_original_part_a2.py` copies the user's original Rmd into `analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd` and replaces **only** the original `SEN <- read_csv(file_path, ...)` import to read the corrected snapshot. Its source checks prove the rest of the notebook is byte-for-byte preserved after undoing that one replacement. SHA-256 of supplied original is recorded.

**This runner intentionally does NOT execute a substituted modern plot pipeline or rewrite analytical chunks.** The generated working copy also normalizes five source-plot variable references from `SEN3` to `SEN`; the SHA-locked 30-chunk original remains intact. The notebook still has original working-directory assumptions, manual intermediary files and exploratory sections. It must be run/migrated with an output-by-output fidelity check before reporting a full reproducible Part A2 run. The 3,434-genome selected cohort must not be silently substituted for the original 33,831 master.
