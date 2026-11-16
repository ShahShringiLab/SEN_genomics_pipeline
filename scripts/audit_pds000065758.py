#!/usr/bin/env python3
"""Independent SEN SNP-cluster/gyrA-D87Y verification (stdlib only).

Inputs are the archived exact repo files; quoted-field normalization mirrors
workflow/A2_Genotype/lib/common.R without changing the source CSV on disk.
This program deliberately prints raw and cross-checked counts, not presumed values.
"""
from __future__ import annotations

import csv
import io
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = "PDS000065758"
PLACEHOLDER = "__SEN_LITERAL_COMMA__"
D87Y = re.compile(r"(?<![A-Za-z0-9])gyrA[_-]?D87Y(?![A-Za-z0-9])", re.I)

def normalized_rows(path: Path, expected_rows: int, expected_cols: int):
    text = path.read_text(encoding="utf-8-sig")
    assert PLACEHOLDER not in text
    escape_count = text.count(r"\,")
    malformed_exudate = text.count(',""Exudate, Left Heel"","clinical",')
    text = text.replace(r"\,", PLACEHOLDER)
    text = text.replace(
        ',""Exudate, Left Heel"","clinical",',
        ',"Exudate, Left Heel","clinical",',
    )
    reader = csv.DictReader(io.StringIO(text, newline=""))
    header = reader.fieldnames or []
    assert len(header) == expected_cols, (path, len(header), header)
    records = []
    anomalous = []
    for i, row in enumerate(reader, start=2):
        if None in row or any(v is None for v in row.values()):
            anomalous.append(i)
        for key, value in row.items():
            if value is not None:
                row[key] = value.replace(PLACEHOLDER, ",")
        records.append(row)
    assert len(records) == expected_rows, (path, len(records), expected_rows)
    assert not anomalous, ("Incorrect CSV field alignment", path, anomalous[:20])
    print(f"[PASS] {path.name}: {len(records):,} rows x {len(header)} cols; "
          f"normalized {escape_count} escaped commas; {malformed_exudate} doubled-quote case(s)")
    return records

def cluster(row):
    s = str(row.get("SNP cluster", "") or "").strip()
    return re.sub(r"\.\d+$", "", s).upper()

def is_d87y(row, col):
    return bool(D87Y.search(str(row.get(col, "") or "")))

def tabulate(rows, label):
    print(f"\n=== {label} ===")
    run_ids = [r.get("Run") for r in rows]
    print(f"Total genomes: {len(rows):,}")
    print(f"Unique runs: {len(set(run_ids)):,}")
    assert len(set(run_ids)) == len(rows)
    members = [r for r in rows if cluster(r) == TARGET]
    print(f"Cluster {TARGET}: {len(members):,}")
    for col in ("AMR genotypes core", "AMR genotypes"):
        if col not in rows[0]:
            print(f"[INFO] Column absent: {col}")
            continue
        positives = [r for r in rows if is_d87y(r,col)]
        both = [r for r in members if is_d87y(r,col)]
        print(f"{col} | D87Y-positive = {len(positives):,}; "
              f"in {TARGET} = {len(both):,}; "
              f"fraction = {100*len(both)/len(positives):.4f}%; "
              f"outside = {len(positives)-len(both):,}")
    return members

def main():
    master = normalized_rows(
        ROOT/"data/part_a/genotype/US-SEN_12192025.csv", 33831, 51
    )
    clades = {}
    with (ROOT/"metadata/final_clade_metadata.tsv").open(newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            clades[r["Sample"]] = r["Lineage"]
    assert Counter(clades.values())["Clade 3"] == 1357
    print(f"[PASS] Clade assignments: {len(clades)} including reference")
    members_full = tabulate(master, "Full 33,831-genome NCBI metadata")

    selected = normalized_rows(
        ROOT/"data/part_a/selection/Finalupdated_SEN_Genomes1.csv", 3434, 69
    )
    selected = [r for r in selected if r.get("Run") in clades and
                clades.get(r.get("Run")) != "Reference"]
    assert len(selected) == 3306, len(selected)
    members_selected = tabulate(selected, "Final selected 3,306-genome phylogeny")
    counts = Counter(clades.get(r["Run"]) for r in members_selected)
    clade3 = [r for r in selected if clades.get(r["Run"]) == "Clade 3"]
    assert len(clade3) == 1357, len(clade3)
    print(f"Cluster membership by clade: {dict(counts)}")
    print(f"Clade 3 in {TARGET}: {counts.get('Clade 3',0)}/1357 = "
          f"{100*counts.get('Clade 3',0)/1357:.6f}%")
    print(f"Cluster members outside Clade 3: "
          f"{sum(n for c,n in counts.items() if c != 'Clade 3')}")
    print(f"Clade 3 members in other/missing clusters: "
          f"{1357-counts.get('Clade 3',0)}")
    unexpected = [r["Run"] for r in clade3 if cluster(r) != TARGET]
    print(f"Off-cluster Clade 3 SRRs ({len(unexpected)}): {', '.join(unexpected)}")

    # All selected records must agree with the full master on SNP-cluster membership.
    selected_full = {r["Run"]:r for r in master if r["Run"] in clades}
    assert len(selected_full) == 3306
    mismatches = [(r["Run"], cluster(r), cluster(selected_full[r["Run"]]))
                  for r in selected if cluster(r) != cluster(selected_full[r["Run"]])]
    print(f"SNP-cluster inconsistencies between corrected selected and full master: {len(mismatches)}")
    if mismatches:
        print(mismatches[:25])

if __name__ == "__main__":
    main()
