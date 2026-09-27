#!/usr/bin/env python3
import os
import re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
meta_path = Path(os.environ.get("SEN_METADATA_FILE", ROOT / "metadata" / "SEN_Genomes.csv"))
out_path = Path(os.environ.get("SEN_SRR_LIST", ROOT / "SRR_list.clean.txt"))

df = pd.read_csv(meta_path, low_memory=False, encoding="utf-8-sig")
if "Run" not in df.columns:
    raise SystemExit("[ERROR] metadata file lacks Run column")

runs = []
for value in df["Run"]:
    m = re.search(r"SRR\d+", str(value))
    if m:
        runs.append(m.group(0))

runs = sorted(set(runs))
if len(runs) != 3434:
    raise SystemExit(f"[ERROR] expected 3434 unique SRRs, found {len(runs)}")

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text("\n".join(runs) + "\n")
print(f"[PASS] wrote {len(runs)} SRRs -> {out_path}")
