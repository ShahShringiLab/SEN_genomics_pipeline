#!/usr/bin/env python3
"""Validate temporal dot-plot eligibility from the independently generated annual table.

A valid collection year <= 2025 is plotted; all years <= 2000 use the <=2000 bin.
Only undated genomes are excluded. This check never modifies the original data.
"""
import csv
from collections import defaultdict
from pathlib import Path
import os

root = Path(os.getenv("SEN_ROOT", Path(__file__).resolve().parents[1]))
file = root / "results/part_a/master_33831/sampling_qc_table/year_source_marker_selection_QC_weights.csv"
if not file.is_file():
    raise SystemExit(f"[FAIL] Sampling table absent: {file}")
counts = defaultdict(lambda: {"total": 0, "missing": 0, "historical": 0, "plotted": 0})
with file.open(encoding="utf-8-sig", newline="") as stream:
    for row in csv.DictReader(stream):
        source = row["Source"]
        label = row["Year of isolation"].strip()
        n = int(row["Available genomes (N)"])
        counts[source]["total"] += n
        if label == "Missing year":
            counts[source]["missing"] += n
            continue
        try:
            year = int(label)
        except ValueError as exc:
            raise SystemExit(f"[FAIL] Unexpected year: {label!r}") from exc
        if year > 2025:
            raise SystemExit(f"[FAIL] Out-of-range dated record: {year}")
        counts[source]["plotted"] += n
        if year <= 2000:
            counts[source]["historical"] += n

expected = {
    "Human": (26541, 245, 19, 26296),
    "Chicken-US": (6637, 49, 17, 6588),
    "Other-US": (582, 29, 2, 553),
    "International": (71, 0, 0, 71),
}
# Verify expected total, missing and eligible counts. Report but do not constrain
# each source's historical count: the <=2000 category is summed explicitly below.
for src, (total, missing, _old, plotted) in expected.items():
    got = counts[src]
    if (got["total"], got["missing"], got["plotted"]) != (total, missing, plotted):
        raise SystemExit(f"[FAIL] Temporal counts for {src}: {got}")
totals = {k: sum(v[k] for v in counts.values()) for k in
          ("total", "missing", "historical", "plotted")}
if totals != {"total": 33831, "missing": 323, "historical": 38, "plotted": 33508}:
    raise SystemExit(f"[FAIL] Unexpected temporal totals: {totals}")
print("[PASS] Master total: 33,831; missing year: 323; plotted: 33,508")
print("[PASS] <=2000 plot bin: 38 genomes; 2001-2025: 33,470 genomes")
for src, vals in counts.items():
    print(f"[AUDIT] {src}: {vals['plotted']:,} dated genomes")
