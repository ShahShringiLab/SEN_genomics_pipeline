#!/usr/bin/env python3
"""Merge verified QC edits into complete 3,434-isolate SEN metadata.

Non-destructive: original 69-column metadata and 11-column QC audit are never
overwritten. Run+Isolate is the exact join key; output retains every isolate,
all original metadata columns, original QR calls, corrected QR/FQR calls,
clade assignments, QC inclusion status, and per-sample audit/rejection notes.
"""
from __future__ import annotations

import csv
import json
import os
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("SEN_ROOT", Path(__file__).resolve().parents[1]))
CANONICAL = ROOT / "data/part_a/selection/Finalupdated_SEN_Genomes1.csv"
QC_FILE = ROOT / "data/part_a/selection/Selected_US-SEN_12192025_QCAUDIT_R1_qr_clade.csv"
CLADE_FILE = ROOT / "metadata/final_clade_metadata.tsv"
OUT = ROOT / "results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv"
SUMMARY = ROOT / "results/part_a/cohort/QC_merged_3434_validation.json"


def load(path: Path, delimiter: str = ",") -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        assert reader.fieldnames, f"Missing header: {path}"
        rows = list(reader)
        assert all(None not in row for row in rows), f"Malformed CSV: {path}"
        return reader.fieldnames, rows


def key(row: dict[str, str]) -> tuple[str, str]:
    return row["Run"].strip(), row["Isolate"].strip()


meta_cols, metadata = load(CANONICAL)
qc_cols, audit = load(QC_FILE)
_, clade_rows = load(CLADE_FILE, "\t")
assert len(metadata) == len(audit) == 3434
assert len(meta_cols) == 69, f"Expected complete 69-column original metadata, got {len(meta_cols)}"
assert len({key(x) for x in metadata}) == len({key(x) for x in audit}) == 3434
assert {key(x) for x in metadata} == {key(x) for x in audit}
assert len({row["Run"] for row in metadata}) == 3434
assert len({row["Isolate"] for row in metadata}) == 3434
assert len({row["Run"] for row in audit}) == 3434
assert set(qc_cols) == {"Run", "Isolate", "QR_genotype_original", "QR_genotype",
                       "QR_status", "QR_marker_count", "QR_markers",
                       "QR_corrected", "Clade", "Final_phylogeny", "AUDIT"}
by_key = {key(row): row for row in audit}
clade_map = {row["Sample"]: row["Lineage"] for row in clade_rows
             if row["Sample"] != "Reference" and row["Lineage"] not in ("", "NA")}
assert len(clade_map) == 3306

merged = []
corrected = []
for original in metadata:
    q = by_key[key(original)]
    assert original["QR_genotype"] == q["QR_genotype_original"], (
        f"Archived QR call differs for {key(original)}"
    )
    old = int(original["QR_genotype"])
    new = int(q["QR_genotype"])
    assert old in (0, 1) and new in (0, 1)
    assert q["QR_corrected"] == ("YES" if old != new else "NO")
    assert q["QR_status"] == ("QR/FQR marker-positive" if new else "QR/FQR marker-negative")
    assert (int(q["QR_marker_count"]) > 0) == bool(new)
    assert q["Final_phylogeny"] in ("Included", "Excluded")
    expected_clade = clade_map.get(q["Run"], "")
    assert q["Clade"] == expected_clade, (
        f"QC and final-clade metadata disagree for {key(original)}: "
        f"{q['Clade']!r} vs {expected_clade!r}"
    )
    assert (q["Final_phylogeny"] == "Included") == bool(expected_clade)
    if new != old:
        corrected.append({"Run": q["Run"], "Isolate": q["Isolate"],
                          "from": old, "to": new})
    row = dict(original)
    row["QR_genotype_archived"] = original["QR_genotype"]
    row["FQR_genotype_archived"] = original["FQR_genotype"]
    row["QR_genotype"] = q["QR_genotype"]
    row["FQR_genotype"] = q["QR_genotype"]
    for column in ("QR_status", "QR_marker_count", "QR_markers",
                   "QR_corrected", "Clade", "Final_phylogeny", "AUDIT"):
        row[column] = q[column]
    merged.append(row)

assert len(corrected) == 3, f"Expected 3 QR edits, found {len(corrected)}"
assert {x["Isolate"] for x in corrected} == {
    "PDT001102834.1", "PDT001168391.1", "PDT001168656.1"
}
counts = Counter(row["Final_phylogeny"] for row in merged)
assert counts == Counter({"Included": 3306, "Excluded": 128})
assert sum(row["AUDIT"].startswith("PASS") for row in merged) == 3306
assert sum(row["AUDIT"].startswith("FAIL") for row in merged) == 128
assert sum(row["AUDIT"].startswith("PASS: Historical coverage exception") for row in merged) == 2
assert sum(row["QR_genotype"] == "1" for row in merged) == 2288
fields = meta_cols + [
    "QR_genotype_archived", "FQR_genotype_archived", "QR_status",
    "QR_marker_count", "QR_markers", "QR_corrected",
    "Clade", "Final_phylogeny", "AUDIT",
]
assert len(fields) == len(set(fields)), "QC merge would overwrite a preexisting field"
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(merged)

summary = {
    "available_master": 33831,
    "full_selected_metadata": len(metadata),
    "columns_original": len(meta_cols),
    "columns_merged": len(fields),
    "unique_runs": len({row["Run"] for row in merged}),
    "unique_isolates": len({row["Isolate"] for row in merged}),
    "qc_included": counts["Included"],
    "qc_excluded": counts["Excluded"],
    "documented_historical_coverage_exceptions_retained": 2,
    "qr_marker_positive": sum(row["QR_genotype"] == "1" for row in merged),
    "qr_corrections": corrected,
    "metadata_source": str(CANONICAL.relative_to(ROOT)),
    "qc_source": str(QC_FILE.relative_to(ROOT)),
    "clade_reference": str(CLADE_FILE.relative_to(ROOT)),
    "output": str(OUT.relative_to(ROOT)),
}
SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
print("[PASS] QC merged into all 3,434 complete selected metadata records.")
print("[PASS] 3,306 included; 128 excluded; 3 QR/FQR calls corrected.")
print(f"[OUTPUT] {OUT}")
print(f"[OUTPUT] {SUMMARY}")
