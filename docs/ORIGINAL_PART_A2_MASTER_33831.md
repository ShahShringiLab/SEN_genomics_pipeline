# Original Part A2: faithful Rmd preservation (33,831 genomes)

**Authoritative source:** user's supplied `SEN_Genomes2.Rmd`. It is longer than the archived `provenance/SEN_Genomes2.Rmd`; they are **not interchangeable**. The full uploaded original is not committed on this branch. Copy it to the local repository before running.

```bash
git clone --branch feature/part-a2-master-33831-parsed-snapshot-20261021 \
 https://github.com/ShahShringiLab/SEN_genomics_pipeline.git SEN_MasterA2
cd SEN_MasterA2
mkdir -p provenance
cp /path/to/your/uploaded/SEN_Genomes2.Rmd provenance/SEN_Genomes2_uploaded_original_20261021.Rmd
source ~/miniconda3/etc/profile.d/conda.sh
conda install -n sen_stats -c conda-forge r-writexl -y
bash scripts/run_master_part_a2.sh
```

The runner first normalizes CSV quoting/escaped commas and exports **exactly 33,831 rows × original 51 columns** as immutable input snapshot at
`results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.xlsx`
and equivalent CSV, with field-alignment audit and checksums.

Then `scripts/prepare_original_part_a2.py` copies the user's original Rmd into `analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd` and replaces **only** the original `SEN <- read_csv(file_path, ...)` import to read the corrected snapshot. Its source checks prove the rest of the notebook is byte-for-byte preserved after undoing that one replacement. SHA-256 of supplied original is recorded.

**This runner intentionally does NOT execute a substituted modern plot pipeline or rewrite analytical chunks.** The notebook still has original working-directory assumptions, manual intermediary files and exploratory sections. It must be run/migrated with an output-by-output fidelity check before reporting a full reproducible Part A2 run. The 3,434-genome selected cohort must not be silently substituted for the original 33,831 master.
