#!/usr/bin/env python3
import csv
import os
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(os.environ.get("SEN_SMOKE_ROOT", ROOT / "tests" / "smoke_test" / "work"))
MANIFEST = ROOT / "tests" / "smoke_test" / "smoke_manifest.tsv"
FULL_META = ROOT / "metadata" / "SEN_Genomes.csv"

SMOKE_META_DIR = WORK / "metadata"
ITOL_DIR = WORK / "ITOL"
SNIPPY_DIR = WORK / "Snippy_output"
SNP_MASTER = SNIPPY_DIR / "SENBIO_RECOVERED_REPORT.csv"
SOURCE_ITOL = WORK / "itol_1_source.txt"
YEAR_ITOL = WORK / "itol_4_collection_year.txt"

SMOKE_META_DIR.mkdir(parents=True, exist_ok=True)
ITOL_DIR.mkdir(parents=True, exist_ok=True)
SNIPPY_DIR.mkdir(parents=True, exist_ok=True)

manifest = pd.read_csv(MANIFEST, sep="\t")
isolates = manifest.loc[manifest["role"] == "isolate"].copy()

# Canonical smoke clade metadata.
clade = isolates[["sample_id", "clade"]].rename(
    columns={"sample_id": "Sample", "clade": "Lineage"}
)
clade.to_csv(SMOKE_META_DIR / "final_clade_metadata.tsv", sep="\t", index=False)

# Legacy clade-list fixtures used only to exercise 14_clademetadata.py.
clade1 = isolates.loc[isolates["clade"].isin(["Clade 1A", "Clade 1B"]), "sample_id"]
clade1b = isolates.loc[isolates["clade"] == "Clade 1B", "sample_id"]
clade2 = isolates.loc[isolates["clade"] == "Clade 2", "sample_id"]
clade3 = isolates.loc[isolates["clade"] == "Clade 3", "sample_id"]

for name, series in [
    ("Clade1.txt", clade1),
    ("Clade1b.txt", clade1b),
    ("Clade2.txt", clade2),
    ("Clade3.txt", clade3),
]:
    (ITOL_DIR / name).write_text("\n".join(series.astype(str)) + "\n")

# Source iTOL fixture from the smoke manifest.
with SOURCE_ITOL.open("w") as fh:
    fh.write("DATASET_COLORSTRIP\n")
    fh.write("SEPARATOR COMMA\n")
    fh.write("DATASET_LABEL,Source\n")
    fh.write("COLOR,#000000\n")
    fh.write("DATA\n")
    for _, row in isolates.iterrows():
        fh.write(f"{row['sample_id']},#000000,{row['source_group']}\n")

# Collection-year fixture. Prefer a year from canonical metadata if one can be
# recovered robustly; otherwise use Unknown, which the analysis supports.
meta = pd.read_csv(FULL_META, low_memory=False, encoding="utf-8-sig")
run_col = next((c for c in ["Run", "Sample", "Accession", "SRR"] if c in meta.columns), None)
if run_col is None:
    raise SystemExit("[ERROR] Cannot locate run-accession column in SEN_Genomes.csv")

meta["_run"] = meta[run_col].astype(str).str.strip()
meta = meta[meta["_run"].isin(set(isolates["sample_id"]))].copy()

def infer_year(row):
    candidates = [
        c for c in meta.columns
        if re.search(r"(collection.*year|year.*collection|collection.*date|date.*collection|^year$)", c, re.I)
    ]
    for c in candidates:
        val = str(row.get(c, "")).strip()
        m = re.search(r"(19|20)\d{2}", val)
        if m:
            y = int(m.group(0))
            if y <= 2000:
                return "<=2000"
            if y <= 2005:
                return "2001-2005"
            if y <= 2010:
                return "2006-2010"
            if y <= 2015:
                return "2011-2015"
            if y <= 2020:
                return "2016-2020"
            if y <= 2025:
                return "2021-2025"
    return "Unknown"

years = {r["_run"]: infer_year(r) for _, r in meta.iterrows()}
with YEAR_ITOL.open("w") as fh:
    fh.write("DATASET_COLORSTRIP\n")
    fh.write("SEPARATOR COMMA\n")
    fh.write("DATASET_LABEL,Collection year\n")
    fh.write("COLOR,#000000\n")
    fh.write("DATA\n")
    for sample in isolates["sample_id"]:
        fh.write(f"{sample},#000000,{years.get(sample, 'Unknown')}\n")

# Build an analysis-compatible SNP master from per-sample Snippy tables.
frames = []
for sample in isolates["sample_id"]:
    tab = SNIPPY_DIR / sample / "snps.tab"
    if not tab.exists() or tab.stat().st_size == 0:
        raise SystemExit(f"[ERROR] Missing Snippy tab for {sample}: {tab}")

    df = pd.read_csv(tab, sep="\t", low_memory=False)
    if "POS" not in df.columns:
        raise SystemExit(f"[ERROR] Snippy tab lacks POS column: {tab}")
    df.insert(0, "Sample", sample)

    # Normalize common Snippy annotation column names to downstream names.
    rename = {}
    if "GENE" not in df.columns and "GENE" in [str(c).upper() for c in df.columns]:
        for c in df.columns:
            if str(c).upper() == "GENE":
                rename[c] = "GENE"
    if rename:
        df = df.rename(columns=rename)

    frames.append(df)

master = pd.concat(frames, ignore_index=True, sort=False)
for col in ["REF", "ALT", "GENE", "EFFECT", "AA_CHANGE"]:
    if col not in master.columns:
        master[col] = pd.NA
master.to_csv(SNP_MASTER, index=False)

print(f"[PASS] smoke clade metadata: {SMOKE_META_DIR / 'final_clade_metadata.tsv'}")
print(f"[PASS] smoke source iTOL:    {SOURCE_ITOL}")
print(f"[PASS] smoke year iTOL:      {YEAR_ITOL}")
print(f"[PASS] smoke SNP master:     {SNP_MASTER} ({len(master)} rows)")
