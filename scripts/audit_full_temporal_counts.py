#!/usr/bin/env python3
"""Validate the 1990-2025 sampling audit and reconcile it with full-master plot totals.

The sampling/QC table deliberately excludes 18 pre-1990 records; full-master
plots retain these genomes in the <=2000 bin. Do not equate the denominators.
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
    "Human": (26525, 245, 9, 26280),
    "Chicken-US": (6637, 49, 1, 6588),
    "Other-US": (580, 29, 9, 551),
    "International": (71, 0, 1, 71),
}
# Verify expected total, missing and eligible counts. Report but do not constrain
# each source's historical count: the <=2000 category is summed explicitly below.
for src, (total, missing, historical, plotted) in expected.items():
    got = counts[src]
    if (got["total"], got["missing"], got["historical"], got["plotted"]) != (total, missing, historical, plotted):
        raise SystemExit(f"[FAIL] Temporal counts for {src}: {got}")
totals = {k: sum(v[k] for v in counts.values()) for k in
          ("total", "missing", "historical", "plotted")}
if totals != {"total": 33813, "missing": 323, "historical": 20, "plotted": 33490}:
    raise SystemExit(f"[FAIL] Unexpected temporal totals: {totals}")
print("[PASS] Sampling audit (1990-2025 plus undated): 33,813; missing year: 323; eligible dated: 33,490")
print("[PASS] Sampling-table <=2000 bin: 20 genomes; 2001-2025: 33,470")
# Separate documented master-population plot scope: add the 18 pre-1990
# genomes excluded from the sampling table (Human=16, Other-US=2).
pre_1990_by_source = {"Human": 16, "Chicken-US": 0, "Other-US": 2, "International": 0}
full_master_total = totals["total"] + sum(pre_1990_by_source.values())
full_master_dated = totals["plotted"] + sum(pre_1990_by_source.values())
full_master_historical = totals["historical"] + sum(pre_1990_by_source.values())
if (full_master_total, full_master_dated, full_master_historical) != (33831, 33508, 38):
    raise SystemExit("[FAIL] Sampling-to-master denominator reconciliation failed")
print("[PASS] Full-master plot reconciliation: 33,831 total; 33,508 dated; <=2000 bin 38")
for src, vals in counts.items():
    print(f"[AUDIT] {src}: {vals['plotted']:,} dated genomes")
