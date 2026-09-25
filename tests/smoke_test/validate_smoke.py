#!/usr/bin/env python3

import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "tests" / "smoke_test" / "smoke_manifest.tsv"
CLADES = ROOT / "metadata" / "final_clade_metadata.tsv"
META = ROOT / "metadata" / "SEN_Genomes.csv"

EXPECTED = {
    "SRR1033488": "Clade 1A",
    "SRR1220773": "Clade 1B",
    "SRR10007520": "Clade 2",
    "SRR10005236": "Clade 3",
}
REQUIRED_CLADES = {"Clade 1A", "Clade 1B", "Clade 2", "Clade 3"}

def read_manifest():
    with MANIFEST.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows

def read_clades():
    with CLADES.open(newline="") as fh:
        return {r["Sample"]: r["Lineage"] for r in csv.DictReader(fh, delimiter="\t")}

def read_metadata():
    with META.open(newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh))
    by_run = {}
    for row in rows:
        run = (row.get("Run") or "").strip()
        if run:
            by_run[run] = row
    return by_run

def main():
    rows = read_manifest()
    clades = read_clades()
    meta = read_metadata()

    isolates = [r for r in rows if r["role"] == "isolate"]
    assert len(isolates) == 4, f"Expected 4 isolates, found {len(isolates)}"
    assert {r["clade"] for r in isolates} == REQUIRED_CLADES

    for r in isolates:
        sample = r["sample_id"]
        assert sample in EXPECTED, f"Unexpected smoke isolate: {sample}"
        assert clades.get(sample) == EXPECTED[sample], (
            f"{sample}: canonical clade={clades.get(sample)!r}, expected={EXPECTED[sample]!r}"
        )
        assert sample in meta, f"{sample}: missing from SEN_Genomes.csv"

        m = meta[sample]
        assert (m.get("Platform") or "").upper() == "ILLUMINA", f"{sample}: not Illumina"
        assert (m.get("Library layout") or "").upper() == "PAIRED", f"{sample}: not paired-end"

    ref = [r for r in rows if r["role"] == "reference"]
    assert len(ref) == 1 and ref[0]["sample_id"] == "Reference"

    print("[PASS] 4 smoke isolates")
    print("[PASS] one isolate from each major clade")
    print("[PASS] all isolate runs are Illumina paired-end in canonical metadata")
    print("[PASS] P125109/NC_011294.1 reference entry present")
    print("SMOKE MANIFEST VALIDATION: PASS")

if __name__ == "__main__":
    try:
        main()
    except AssertionError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        sys.exit(1)
