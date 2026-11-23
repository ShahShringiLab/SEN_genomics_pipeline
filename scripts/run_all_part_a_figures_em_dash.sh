#!/usr/bin/env bash
# Regenerate the complete original Part A figure set, then collect TIFFs only.
# This runs the existing validated full Part A workflow, including all
# original full-master and selected-cohort analysis chunks and add-on figures.
# Writes ONLY within this checkout. Does not alter external snapshot folders.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"

# Support non-initialized Conda shells, including AjuloPC WSL.
if ! command -v conda >/dev/null 2>&1; then
  for candidate in "$HOME/miniconda3/bin/conda" "$HOME/miniforge3/bin/conda" \
                   "$HOME/anaconda3/bin/conda"; do
    if [[ -x "$candidate" ]]; then
      export PATH="$(dirname "$candidate"):$PATH"
      break
    fi
  done
fi
command -v conda >/dev/null 2>&1 || { echo "[ERROR] Conda not found" >&2; exit 1; }

echo "[INFO] Running complete native Part A plotting workflow (not a 7-figure subset)."
bash scripts/run_part_a.sh

echo "[INFO] Rebuilding manuscript Figure 2, 3, S2, S3 with focused label checks."
SEN_PARTA_SNAPSHOTS="$ROOT" bash scripts/run_part_a_label_figures.sh

python3 - <<'PY'
from pathlib import Path
import shutil

root = Path.cwd()
base = root / "results" / "part_a"
out = base / "all_figures_em_dash_tiff"
out.mkdir(parents=True, exist_ok=True)
sources = sorted(p for p in base.rglob("*") if p.is_file()
                 and p.suffix.lower() in (".tif", ".tiff")
                 and out not in p.parents)
if not sources:
    raise SystemExit("[ERROR] No Part A TIFFs were generated")
manifest = []
for src in sources:
    relative = src.relative_to(base)
    dest = out / relative
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    manifest.append(str(relative))
required = [
    "FQR_marker_counts_top100",
    "FQ_burden",
    "FQ_burden_5yr",
    "FQ_QRFQR_vs_QFQsusceptible_PERCENT_by_year",
    "FQ_QRFQR_vs_QFQsusceptible_COUNT_by_year",
    "FQ_QSFQS_TOP5_QRFQR_PERCENT_by_year",
    "FQ_QSFQS_TOP5_QRFQR_COUNT_by_year",
]
corrected = out / "mdpi_label_corrections"
missing = [name for name in required if not (corrected / (name + ".tif")).is_file()]
if missing:
    raise SystemExit("[ERROR] Missing mandatory corrected TIFFs: " + ", ".join(missing))
(out / "TIFF_MANIFEST.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
print(f"[PASS] {len(sources)} Part A TIFF exports collected separately (including all 7 mandatory corrected manuscript plots).")
print(f"[OUTPUT] {out}")
PY
