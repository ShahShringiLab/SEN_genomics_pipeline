#!/usr/bin/env python3
"""Fail-closed independent audit of original full-master results and selected QC.

Never substitute a 3,434-genome denominator for master-wide analyses.
Do not silently change reported source counts to match expectations.
"""
from pathlib import Path
from collections import Counter
import csv
import json
import os

ROOT = Path(os.environ.get("SEN_ROOT", Path(__file__).resolve().parents[1]))
FULL = ROOT / "results/part_a/full_master_original"
SELECTED = ROOT / "results/part_a/selected_original"
COHORT = ROOT / "results/part_a/cohort"


def read_csv(name):
    p = Path(name)
    if not p.is_file() or p.stat().st_size == 0:
        raise SystemExit(f"[FAIL] Required original result absent: {p}")
    with p.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise SystemExit(f"[FAIL] Original result is empty: {p}")
    return rows


def require(condition, message):
    if not condition:
        raise SystemExit(f"[FAIL] {message}")


for filename in (
    "FQR_marker_counts_top100.png",
    "FQR_marker_counts_top100.tif",
    "FQ_burden.png",
    "FQ_burden_5yr.png",
    "QR_diversity_and_composition_side_by_side.png",
):
    require((FULL / filename).is_file() and (FULL / filename).stat().st_size > 0,
            f"Full-master original figure missing: {filename}")

data = read_csv(FULL / "SEN_TableS1_final_source_decomposition.csv")
expected = {
    ("International", "Others"): 67,
    ("International", "Chicken"): 4,
    ("Other-US", "Others"): 448,
    ("Other-US", "Cattle"): 60,
    ("Other-US", "Turkey"): 48,
    ("Other-US", "Swine"): 21,
    ("Other-US", "Sheep/goat"): 5,
}
observed = {}
for row in data:
    key = (row["Source category"].strip(), row["Decomposed source information"].strip())
    require(key not in observed, f"Duplicate source classification row: {key}")
    observed[key] = int(float(row["Number of genomes (n)"]))
require(observed == expected,
        "Correctly parsed full-master source distribution differs from original Table S1. "
        f"Expected {expected}, observed {observed}. DO NOT silently overwrite the table.")
require(sum(n for (g, _), n in observed.items() if g == "International") == 71,
        "International source total differs from 71")
require(sum(n for (g, _), n in observed.items() if g == "Other-US") == 582,
        "Other-US source total differs from 582")
print("[PASS] ORIGINAL FULL MASTER Table S1 matches 71 International and 582 Other-US genomes")
print("[PASS] International Others=67 Chicken=4; Other-US Others=448 Cattle=60 Turkey=48 Swine=21 Sheep/goat=5")

markers = read_csv(FULL / "FQR_marker_frequency_percentage_95CI_with_susceptible.csv")
susceptible = [r for r in markers if r["token_pretty"] == "Q/FQ-susceptible"]
require(len(susceptible) == 1, "Missing unique no-marker category from full master")
n_no_marker = int(float(susceptible[0]["n_isolates"]))
print(f"[AUDIT] Original full-master QR/FQR markers: {33831-n_no_marker:,} positive; "
      f"{n_no_marker:,} without detected markers; N=33,831")
require(n_no_marker + (33831-n_no_marker) == 33831, "Master marker denominator changed")

merged = read_csv(COHORT / "Finalupdated_SEN_Genomes1_QC_merged.csv")
require(len(merged) == 3434, f"Selected QC metadata has {len(merged)} rows, expected 3434")
require(len(merged[0]) == 78, "QC metadata must preserve 69 originals + 9 QC fields")
require(len({r["Run"] for r in merged}) == 3434, "Selected QC Run values are not unique")
require(len({r["Isolate"] for r in merged}) == 3434, "Selected QC Isolate values are not unique")
require(Counter(r["Final_phylogeny"] for r in merged) ==
        Counter({"Included": 3306, "Excluded": 128}),
        "Selected QC inclusion does not reconcile to 3306/128")
require(sum(r["QR_corrected"] == "YES" for r in merged) == 3, "QR correction count != 3")
require(sum(r["AUDIT"].startswith("PASS: Historical coverage exception") for r in merged) == 2,
        "Expected 2 documented QC historical-coverage exceptions")
print("[PASS] Selected 3,434 original metadata records retained, QC 3,306/128, 3 QR corrections")

# The original selected-period figures must not overwrite full-master output.
require((SELECTED / "SEN_Total_TemporalBlocks_QRFQR_Count.png").is_file(),
        "Original selected-versus-total temporal-block figure missing")

ROOT.joinpath("results/part_a/full_master_original/REPRO_AUDIT.json").write_text(
    json.dumps({
        "master_genomes": 33831,
        "dated_for_temporal_plot": 33508,
        "excluded_from_temporal_plot_missing_year": 323,
        "source_table": [{
            "source": k[0], "animal_source": k[1], "n": v
        } for k, v in sorted(observed.items())],
        "marker_positive": 33831-n_no_marker,
        "without_detected_marker": n_no_marker,
        "selected": 3434,
        "qc_included": 3306,
        "qc_excluded": 128,
        "qr_corrections": 3
    }, indent=2)+"\n", encoding="utf-8")
print("[PASS] Original 33,831-master analytical outputs and 3,434-cohort QC fully separated")
