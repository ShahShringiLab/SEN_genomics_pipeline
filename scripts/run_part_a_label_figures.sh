#!/usr/bin/env bash
# Rerender ONLY Figure 2, Figure 3, and the original four separate S2/S3
# panels. Do not modify inputs, original Rmd, or previously saved figures.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"

CONDA_BIN="$(command -v conda 2>/dev/null || true)"
if [[ -z "$CONDA_BIN" ]]; then
  for candidate in "$HOME/miniforge3/bin/conda" "$HOME/miniforge3/condabin/conda" \
                   "$HOME/miniconda3/bin/conda" "$HOME/miniconda3/condabin/conda"; do
    if [[ -x "$candidate" ]]; then CONDA_BIN="$candidate"; break; fi
  done
fi
if [[ -z "$CONDA_BIN" ]]; then
  echo "[ERROR] Cannot find installed Conda/Miniforge." >&2
  exit 1
fi

snapshot_root="${SEN_PARTA_SNAPSHOTS:-$ROOT}"
if [[ ! -d "$snapshot_root" ]]; then
  echo "[ERROR] Snapshot directory not found: $snapshot_root" >&2
  exit 1
fi
snapshot_root="$(cd "$snapshot_root" && pwd)"
export SEN_PARTA_SNAPSHOTS="$snapshot_root"
raw="$snapshot_root/results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv"
derived="$snapshot_root/results/part_a/full_master_original/SEN_33831_original_Rmd_derived.csv"
for file in "$raw" "$derived"; do
  if [[ ! -s "$file" ]]; then
    echo "[ERROR] Required validated original Part A snapshot missing: $file" >&2
    echo "Run original bash scripts/run_part_a.sh once to create the scientific inputs." >&2
    exit 1
  fi
done

# Regenerate the prepared execution copy from the frozen original when absent.
# This does not modify any original manuscript source or scientific data.
if [[ ! -s workflow/A2_Genotype/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd ]]; then
  python3 scripts/prepare_original_part_a2.py
fi

python3 scripts/build_part_a_label_figures.py

out="results/part_a/mdpi_label_corrections"
mkdir -p "$out"
echo "[INFO] Replotting only publication figures using validated snapshots from $snapshot_root."
"$CONDA_BIN" run --no-capture-output -n sen_stats \
  env SEN_ROOT="$ROOT" \
  Rscript "$out/run_figures_2_3_s2_s3.R" 2>&1 | tee "$out/run.log"

echo "[INFO] Exporting seven independent PDFs (retaining any native Figure 2 PDF)."
"$CONDA_BIN" run --no-capture-output -n sen_stats \
  python scripts/export_part_a_figures_pdf.py "$out" 2>&1 | tee -a "$out/run.log"

echo "[PASS] Corrected Figure 2, Figure 3 and four supplementary panels; seven PDFs available."
echo "[OUTPUT] $ROOT/$out/"
