#!/usr/bin/env python3
import os
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sen_root = Path(os.environ.get("SEN_ROOT", ROOT))
meta_path = Path(os.environ.get("SEN_METADATA_FILE", ROOT / "metadata" / "SEN_Genomes.csv"))
source_out = Path(os.environ.get("SEN_SOURCE_ITOL", sen_root / "itol_1_source.txt"))
year_out = Path(os.environ.get("SEN_COLLECTION_YEAR_ITOL", sen_root / "itol_4_collection_year.txt"))

df = pd.read_csv(meta_path, low_memory=False, encoding="utf-8-sig")
for col in ["Run", "Source", "Collection_year"]:
    if col not in df.columns:
        raise SystemExit(f"[ERROR] metadata lacks required column: {col}")

def year_bin(value):
    try:
        y = int(float(value))
    except Exception:
        return "Unknown"
    if y <= 2000: return "<=2000"
    if y <= 2005: return "2001-2005"
    if y <= 2010: return "2006-2010"
    if y <= 2015: return "2011-2015"
    if y <= 2020: return "2016-2020"
    if y <= 2025: return "2021-2025"
    return "Unknown"

def write_colorstrip(path, label, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("DATASET_COLORSTRIP\nSEPARATOR COMMA\n")
        fh.write(f"DATASET_LABEL,{label}\nCOLOR,#000000\nDATA\n")
        for run, value in rows:
            fh.write(f"{run},#000000,{value}\n")

source_rows = []
year_rows = []
seen = set()
for _, r in df.iterrows():
    run = str(r["Run"]).strip()
    if not run.startswith("SRR") or run in seen:
        continue
    seen.add(run)
    source = str(r["Source"]).strip()
    if source.lower() in {"", "nan", "none"}:
        source = "Unknown"
    source_rows.append((run, source))
    year_rows.append((run, year_bin(r["Collection_year"])))

write_colorstrip(source_out, "Source", source_rows)
write_colorstrip(year_out, "Collection year", year_rows)
print(f"[PASS] source covariate: {source_out}")
print(f"[PASS] year covariate:   {year_out}")
