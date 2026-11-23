# Figure 2, Figure 3, S2 and S3: publication typography only

This branch makes precisely the requested presentation changes to the **full 33,831-genome Part A manuscript figures** while retaining the canonical original R Markdown scientific source, selection, all statistics, plot data, palettes, dimensions and axis scales.

## One command (Ubuntu WSL with existing Miniforge)

From the cloned repository (on branch `figure2-figure3-s2-s3-labels-only`):

```bash
bash scripts/run_part_a_label_figures.sh
```

The runner uses the already-generated, validated Part A **raw** and **Rmd-derived** 33,831-row metadata snapshots. If these snapshots are missing, it stops and tells you to run the full Part A workflow; it does not silently replace them or recalculate other analyses. It creates the prepared working Rmd if missing, without altering the source notebook.

If the validated results are in a **different Part A checkout**, you can reuse them without copying the files:

```bash
SEN_PARTA_SNAPSHOTS="$HOME/SEN_PartA_clean" bash scripts/run_part_a_label_figures.sh
```

Replace `~/SEN_PartA_clean` with the path to the checkout that already contains the two validated full-master CSV snapshots. The outputs are still written under the current branch's `results/part_a/mdpi_label_corrections/`.

The run writes to a **new output folder**:
`results/part_a/mdpi_label_corrections/`

## Figure 2 output (same original sizing and colors)

- `FQR_marker_counts_top100.png`
- `FQR_marker_counts_top100.tif`
- `FQR_marker_counts_top100.pdf`

Only displayed axis labels are modified: `Q/FQ-susceptible` becomes `QS/FQS`; gene symbols are italic with mutation suffixes upright; the count text `18375` and `13486` gains commas; the count-axis `n` is italic. Underlying data columns, raw table labels, counts and percentages remain unchanged.

## Figure 3 outputs

- `FQ_burden.png` / `FQ_burden.tif`
- `FQ_burden_5yr.png` / `FQ_burden_5yr.tif` (related 5-year version)

Source facet labels display: Human (*n* = 26,296), Chicken—US (*n* = 6,588), Other—US (*n* = 553), International (*n* = 71). Source names use em dash (U+2014). Gene symbols in marker rows are italic; mutation suffixes are upright. Percent legend range separators become em dash (U+2014); no numerical percent-bin boundary changes.

## Supplementary Figure S2: two **separate** plots

- `FQ_QRFQR_vs_QFQsusceptible_PERCENT_by_year.png` / `.tif`
- `FQ_QRFQR_vs_QFQsusceptible_COUNT_by_year.png` / `.tif`

## Supplementary Figure S3: two **separate** plots

- `FQ_QSFQS_TOP5_QRFQR_PERCENT_by_year.png` / `.tif`
- `FQ_QSFQS_TOP5_QRFQR_COUNT_by_year.png` / `.tif`

The four supplementary panels are exported independently using the original source filenames, **not** assembled into one composite. Their source facet typography is corrected; S3 marker legends also italicize gene symbols. Every original bar position, year, count, percent, color, and plot size is retained.

## Code provenance and change isolation

- Original frozen scientific notebook **unchanged**: `provenance/SEN_Genomes2_uploaded_original_20261021.Rmd`
- Scientific Part A execution generator: `scripts/build_original_a2_master_execution.py`
- Exact, indexed display-only edits: `scripts/part_a_figure_presentation.py`
- Display functions: `scripts/part_a_figure_label_helpers.R`
- Focused figure-only build: `scripts/build_part_a_label_figures.py`
- Runner: `scripts/run_part_a_label_figures.sh`

The original full run (`bash scripts/run_part_a.sh`) also applies the same display-only changes in the generated execution script; original analytical tables and source code are not rewritten. The focused runner avoids rerunning other Part A statistical analyses and never overwrites the previously exported `results/part_a/full_master_original/` figures.

Run status: GitHub source changes are committed, but plots need to be rendered on the user's WSL computer; no local R execution is claimed.


## TIFF-only manuscript figure export

The focused runner retains seven independent TIFF files in
`results/part_a/mdpi_label_corrections/` and removes PNG/PDF companions
from this dedicated directory. No Pillow installation is required.

## Complete Part A figures — em dash and TIFF only

Use the full-workflow runner when **all Part A plots**, not just Figure 2 / Figure 3 / S2 / S3, must be regenerated:

```bash
cd ~/SEN_Figure_R2
git pull origin figure2-figure3-s2-s3-labels-only
bash scripts/run_all_part_a_figures_em_dash.sh
```

This uses the existing frozen source and canonical cohorts, executes the existing Part A
master, selected-genome, temporal Top-5, and International plotting workflows,
and then regenerates the seven specially formatted manuscript figures.
It retains all previous Figure 2/3/S2/S3 changes and converts **en dash (–) to em dash (—) in displayed plot labels** without recoding any underlying category or time period.
Numeric percent-bin cutpoints remain the same.

The consolidated TIFF-only output is:

`results/part_a/all_figures_em_dash_tiff/`

The subdirectories preserve each source figure's original filename and analysis
location. `TIFF_MANIFEST.txt` lists the TIFFs. Some original scripts also create
PNG/PDF files in their normal intermediate folders; **only TIFFs** are collected
in the final folder. The historical scientific Rmd is never altered.

**Important:** The complete runner performs the original Part A computations and
depends on the canonical inputs and packages from `scripts/run_part_a.sh`.
It must not be mistaken for the faster 7-figure-only runner. Re-rendering
and checking its final outputs is performed on the user's WSL machine.
