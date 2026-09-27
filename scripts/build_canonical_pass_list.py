#!/usr/bin/env python3
import os
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
meta_path = Path(os.environ.get("SEN_CLADE_METADATA", ROOT / "metadata" / "final_clade_metadata.tsv"))
out_path = Path(os.environ.get("SEN_PASS_LIST", ROOT / "confirmed_enteritidis_list.txt"))

df = pd.read_csv(meta_path, sep="\t")
required = {"Sample", "Lineage"}
if not required.issubset(df.columns):
    raise SystemExit(f"[ERROR] missing columns: {sorted(required - set(df.columns))}")

study = df[df["Lineage"] != "Reference"].copy()
study = study[study["Sample"].astype(str).str.match(r"^SRR\d+$")].copy()
runs = sorted(study["Sample"].drop_duplicates().tolist())

if len(runs) != 3306:
    raise SystemExit(f"[ERROR] expected 3306 canonical study SRRs, found {len(runs)}")

counts = study["Lineage"].value_counts().to_dict()
expected = {
    "Clade 3": 1357,
    "Clade 2": 894,
    "Clade 1B": 759,
    "Clade 1A": 169,
    "Unassigned": 127,
}
if counts != expected:
    raise SystemExit(f"[ERROR] canonical clade counts differ: {counts}")

out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text("\n".join(runs) + "\n")
print(f"[PASS] wrote canonical 3306-isolate pass list -> {out_path}")
