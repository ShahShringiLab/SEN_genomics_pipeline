#!/usr/bin/env python3
"""Post-run integrity checks for both enriched metadata files and paired stacks."""
import csv
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parents[1]
fullfile = root / "results/part_a/reconciled_metadata/SEN_33831_full_enriched_metadata.csv"
selfile = root / "results/part_a/reconciled_metadata/SEN_3434_selected_enriched_QC_metadata.csv"
countsfile = root / "results/part_a/temporal_top5/SEN_Total_Selected_TemporalBlocks_Top5_QSFQS_data.csv"
rawfile = root / "results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv"
qcfile = root / "results/part_a/cohort/Finalupdated_SEN_Genomes1_QC_merged.csv"

def load(path, expected):
    if not path.is_file() or not path.stat().st_size:
        raise SystemExit(f"[FAIL] Missing enriched output: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        cols, rows = reader.fieldnames, list(reader)
    assert cols and len(rows)==expected,(path, len(rows), expected)
    assert len(rows)==len({(x["Run"],x["Isolate"]) for x in rows}), f"Duplicate Run/Isolate: {path}"
    return cols,rows

raw_cols, raw = load(rawfile,33831)
qc_cols, qc = load(qcfile,3434)
full_cols, full = load(fullfile,33831)
sel_cols, sel = load(selfile,3434)
assert len(raw_cols)==51
assert len(qc_cols)==78
assert all(c in full_cols for c in raw_cols)
assert all(c in sel_cols for c in qc_cols)
for col in ("Rmd_Final_source","Rmd_Primary_source",
            "Audit_Collection_year","Audit_Source","Audit_Temporal_block",
            "Audit_QR_status","Audit_QR_marker_tokens","Audit_QR_marker_count",
            "Audit_Top5_bucket","Audit_Selected_pre_QC",
            "Audit_Analyzed_post_QC"):
    assert col in full_cols and col in sel_cols, col
assert len({r["Run"] for r in full})==33831
assert len({r["Run"] for r in sel})==3434
assert {r["Run"] for r in qc}=={r["Run"] for r in sel}
qc_by_run={r["Run"]:r for r in qc}
selected_by_run={r["Run"]:r for r in sel}
assert sum(r["Audit_Selected_pre_QC"]=="TRUE" for r in full)==3434
assert sum(r["Audit_Analyzed_post_QC"]=="TRUE" for r in full)==3306
assert sum(r["QR_corrected"]=="YES" for r in sel)==3
assert Counter(r["Final_phylogeny"] for r in sel)==Counter({"Included":3306,"Excluded":128})
for r in sel:
    orig=qc_by_run[r["Run"]]
    assert all(r[c]==orig[c] for c in qc_cols), f"Changed selected archival/QC fields: {r['Run']}"
for r in full:
    if r["Run"] in selected_by_run:
        s=selected_by_run[r["Run"]]
        assert r["Audit_Top5_bucket"]==s["Audit_Top5_bucket"]
        assert r["Audit_Source"]==s["Audit_Source"]
        assert r["QC_Final_phylogeny"]==s["Final_phylogeny"]
assert Counter(r["Audit_Source"] for r in full)==Counter(
    {"Human":26541,"Chicken-US":6637,"Other-US":582,"International":71})
valid_cats={"QS/FQS","gyrA_D87Y","gyrA_D87N","gyrA_S83Y","gyrA_D87G","qnrB19","Others"}
assert all(r["Audit_Top5_bucket"] in valid_cats for r in full)
assert sum(bool(r["Audit_Temporal_block"]) for r in full)==33508
assert sum(bool(r["Audit_Temporal_block"]) for r in sel)==3434
chart = root / "results/part_a/temporal_top5"
for population in ("Total_33831","Selected_3434"):
    for unit in ("Count","Percent"):
        stem=f"SEN_{population}_TemporalBlocks_Top5_QSFQS_{unit}"
        for ext in (".png",".tif",".pdf"):
            path=chart/(stem+ext)
            assert path.is_file() and path.stat().st_size>0,path
with countsfile.open(newline="",encoding="utf-8-sig") as fh:
    rows=list(csv.DictReader(fh))
counts=Counter()
percs=Counter()
for r in rows:
    assert r["Category"] in valid_cats
    counts[r["Cohort"]]+=int(r["Count"])
    percs[(r["Cohort"],r["Period"])]+=float(r["Percent"])
assert counts==Counter({"Total available 33,831-genome master":33508,
                       "Selected 3,434-genome cohort":3434}), counts
assert len(percs)==6 and all(abs(p-100)<1e-5 for p in percs.values()),percs
print("[PASS] All original 51 full-master and 78 selected+QC columns are preserved.")
print("[PASS] Enriched master=33,831; selected=3,434; QC included=3,306; excluded=128.")
print("[PASS] Master-to-selected source, QR bucket and QC statuses reconcile by Run.")
print("[PASS] Full temporal stacks=33,508 dated; selected=3,434; exactly 1 marker bucket per genome.")
print("[PASS] Matching count/percentage PNG, TIFF and PDF figures and audit CSVs exported.")
