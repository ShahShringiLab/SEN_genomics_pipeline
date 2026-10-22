#!/usr/bin/env bash
# Rebuild Table S5 safely from sample-tagged AMRFinderPlus master and 3306 clade labels.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export SEN_ROOT="$ROOT"
AMR="${SEN_AMRFINDER_MASTER:-$ROOT/AMRFinderplus/master_AMRFinder_report.tsv}"
CLADE="${SEN_CLADE_METADATA:-$ROOT/metadata/final_clade_metadata.tsv}"
MC="${SEN_MC_RESAMPLES:-100000}"
OUT="$ROOT/iqtree_final/AMRFinder_clade/Split_Clade1A1B/With_Unassigned"

echo "[INPUT] AMRFinder master: $AMR"
echo "[INPUT] Clade metadata:  $CLADE"
echo "[INPUT] Monte Carlo resamples: $MC"
[[ -s "$AMR" ]] || { echo "[FAIL] AMRFinder master missing: $AMR" >&2; exit 1; }
[[ -s "$CLADE" ]] || { echo "[FAIL] Clade metadata missing: $CLADE" >&2; exit 1; }
# Fail before Python analysis if someone supplies the ordinary 22-column
# untagged AMRFinderPlus report. Such rows cannot be linked to phylogenetic tips.
python3 - "$AMR" "$CLADE" <<'PY'
from pathlib import Path
import sys
m=Path(sys.argv[1]); c=Path(sys.argv[2])
with m.open(encoding="utf-8-sig") as h:
    fields=h.readline().rstrip("\r\n").split("\t")
with c.open(encoding="utf-8-sig") as h:
    cf=h.readline().rstrip("\r\n").split("\t")
if not any(x in fields for x in ("Sample","sample","Sample ID","Sample_ID","sample_id")):
    raise SystemExit("[FAIL] Master has no Sample column. Use the sample-tagged master built by workflow/07_amr_vf_plasmid/24_amrfinderplus.sh; never rename Protein id to Sample.")
if not {"Sample","Lineage"}.issubset(cf):
    raise SystemExit("[FAIL] Clade mapping requires Sample and Lineage columns")
print(f"[PASS] AMRFinder columns={len(fields)}; explicit sample ID present; clade metadata schema valid")
PY

echo "[1/2] Test five-clade statistics and significance-letter logic"
conda run --no-capture-output -n sen_stats \
  python tests/test_amrfinder_clade_stats.py

echo "[2/2] Run split clades + Unassigned Table S5, preserving every 0%/100% group"
mkdir -p "$OUT"
conda run --no-capture-output -n sen_stats \
  python analysis/clades/15_amr_byclade.py \
    --amr-master "$AMR" \
    --clade-metadata "$CLADE" \
    --scheme Split_Clade1A1B \
    --mode With_Unassigned \
    --stats-only \
    --mc-resamples "$MC" \
    2>&1 | tee "$OUT/Table_S5_corrected_rerun.log"

echo "[PASS] Corrected Table S5 audit completed."
echo "[RESULTS] $OUT"
echo "[RESULTS] $OUT/AMRFinder_ByClade_Summary.xlsx"
echo "[RESULTS] $OUT/AMRFinder_Omnibus_Stats.tsv"
echo "[RESULTS] $OUT/AMRFinder_Pairwise_Fisher.tsv"
echo "[RESULTS] $OUT/Table_S5_Flagged_Genes_Audit.tsv"
