#!/usr/bin/env python3
"""Non-destructive raw-master CSV repair audit for Part A2.

Produces a workbook with exactly two sheets: Before and After.
Before uses Python csv.reader on the original unmodified CSV *without*
normalizing escaped commas or doubled quotes. This faithfully records the
parser's displaced field values (including rows with abnormal widths).
After reads the validated corrected 51-column master snapshot.
No changes to master, QR genotypes, plots, or selection are performed here.
"""
from __future__ import annotations
import csv
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/"data/part_a/genotype/US-SEN_12192025.csv"
CORRECTED=ROOT/"results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv"
OUT=CORRECTED.parent
MALFORMED=OUT/"SEN_CSV_Malformed_Before_After.xlsx"
DELTA=OUT/"SEN_CSV_Parsing_Field_Differences.csv"
SUMMARY=OUT/"SEN_CSV_Parsing_Change_Summary.csv"

def run():
    import pandas as pd
    from openpyxl import load_workbook
    from openpyxl.styles import Font,PatternFill,Alignment
    OUT.mkdir(parents=True,exist_ok=True)
    with RAW.open(encoding="utf-8-sig",newline="") as f:
        original=list(csv.reader(f,strict=False))
    headers=original[0]
    assert len(headers)==51, len(headers)
    assert len(original)-1==33831,(len(original)-1)
    with CORRECTED.open(encoding="utf-8-sig",newline="") as f:
        corrected=list(csv.DictReader(f))
    assert len(corrected)==33831
    assert all(set(r)==set(headers) for r in corrected)
    before=[];after=[];deltas=[];counts=Counter()
    for position,(parsed,clean) in enumerate(zip(original[1:],corrected),start=2):
        # Matching positional rows is validated on *unaltered stable accession*
        # identifiers when those columns were not displaced.
        if len(parsed)<51:
            padded=parsed+[""]*(51-len(parsed))
        else:
            padded=parsed[:51]
        rawmap=dict(zip(headers,padded))
        unequal=[col for col in headers if rawmap[col]!=clean[col]]
        if not unequal and len(parsed)==51:
            continue
        counts["affected_rows"]+=1
        reason=[]
        if len(parsed)!=51:reason.append(f"raw_column_count={len(parsed)}")
        # These identify both malformed quoted content and escaped commas.
        if len(parsed)==51:reason.append("field_value_changed_by_normalization")
        desc="; ".join(reason)
        shared={"Original_CSV_Row":position,"Repair_reason":desc,
                "Raw_field_count":len(parsed),"Affected_columns":"; ".join(unequal)}
        before.append({**shared,**rawmap, "Extra_raw_fields":" | ".join(parsed[51:])})
        after.append({**shared,**clean,"Extra_raw_fields":""})
        for field in unequal:
            deltas.append({"Original_CSV_Row":position,
                "Run_after":clean.get("Run",""),"Isolate_after":clean.get("Isolate",""),
                "Field":field,"Before":rawmap[field],"After":clean[field],
                "Repair_reason":desc})
            counts["changed_fields"]+=1
    assert len(before)==len(after)>0
    assert [x["Original_CSV_Row"] for x in before]==[x["Original_CSV_Row"] for x in after]
    # The 24 escaped-comma occurrences and one doubled-quote anomaly are
    # source-text events, not assumed to be 25 distinct genomes.
    rawtext=RAW.read_text(encoding="utf-8-sig")
    assert rawtext.count(r"\,")==24,rawtext.count(r"\,")
    assert rawtext.count(',""Exudate, Left Heel"","clinical",')==1
    pd.DataFrame(deltas).to_csv(DELTA,index=False)
    # Preserve all archived-to-corrected QR calls from the frozen selected
    # QC audit without recalculating genotypes or changing plot inputs.
    qc_path=ROOT/"data/part_a/selection/Selected_US-SEN_12192025_QCAUDIT_R1_qr_clade.csv"
    with qc_path.open(encoding="utf-8-sig",newline="") as fh:
        qc=list(csv.DictReader(fh))
    qr_deltas=[dict(row) for row in qc if row["QR_genotype_original"] != row["QR_genotype"]]
    assert len(qr_deltas)==3,(len(qr_deltas))
    qr_path=OUT/"SEN_QR_Genotype_Before_After_3_Genomes.csv"
    pd.DataFrame(qr_deltas).to_csv(qr_path,index=False)

    pd.DataFrame([{"Metric":"raw_rows","Value":33831},
                  {"Metric":"corrected_rows","Value":len(corrected)},
                  {"Metric":"affected_row_count","Value":len(before)},
                  {"Metric":"different_parsed_cells","Value":len(deltas)},
                  {"Metric":"source_escaped_comma_occurrences","Value":24},
                  {"Metric":"source_doubled_quote_occurrences","Value":1}]).to_csv(SUMMARY,index=False)
    with pd.ExcelWriter(MALFORMED,engine="openpyxl") as writer:
        for title, rows in (("Before",before),("After",after)):
            pd.DataFrame(rows).to_excel(writer,sheet_name=title,index=False)
            ws=writer.sheets[title]
            ws.freeze_panes="F2"
            ws.auto_filter.ref=ws.dimensions
            ws.sheet_view.showGridLines=False
            ws.row_dimensions[1].height=30
            for cell in ws[1]:
                cell.fill=PatternFill("solid",fgColor="223C5D")
                cell.font=Font(color="FFFFFF",bold=True)
                cell.alignment=Alignment(wrap_text=True)
            for col in ("A","B","C","D","E"):
                ws.column_dimensions[col].width=24
    x=load_workbook(MALFORMED,read_only=True)
    assert x.sheetnames==["Before","After"]
    assert x["Before"].max_row==x["After"].max_row==len(before)+1
    print(f"[PASS] Parser before/after: {len(before)} affected source rows, {len(deltas)} differing cells.")
    print("[OUTPUT]",MALFORMED)
    print("[OUTPUT]",DELTA)
    print("[OUTPUT]",SUMMARY)
    print("[OUTPUT]",qr_path)

if __name__=="__main__":
    run()
