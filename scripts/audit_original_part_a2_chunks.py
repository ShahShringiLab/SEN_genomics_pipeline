#!/usr/bin/env python3
"""Audit original Part A2 chunk boundaries and optionally instrument a COPY.

No changes to the original notebook, plotting code, or statistical calls.
"""
import argparse,csv,re
from pathlib import Path
root=Path(__file__).resolve().parents[1]
ap=argparse.ArgumentParser()
ap.add_argument("--notebook",type=Path,default=root/"analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd")
ap.add_argument("--instrument",action="store_true",help="Generate a separate guarded notebook copy")
args=ap.parse_args()
lines=args.notebook.read_text(encoding="utf-8").splitlines(keepends=True)
openers=[i for i,x in enumerate(lines) if re.match(r"^\s*```\{r(?:\s|,|\})",x)]
chunks=[]
for number,i in enumerate(openers,1):
    closing=next((j for j in range(i+1,len(lines)) if lines[j].strip()=="```"),None)
    if closing is None:raise SystemExit(f"Unclosed chunk {number} at line {i+1}")
    body="".join(lines[i+1:closing])
    touches=bool(re.search(r"(?m)^\s*SEN\s*(?:<-|=)",body))
    population_risks=[x for x,pat in {
      "row-filter":r"\b(?:filter|slice|distinct|semi_join|inner_join|left_join|right_join|full_join)\s*\(",
      "SEN-reassignment":r"(?m)^\s*SEN\s*(?:<-|=)",
      "disk-input":r"\b(?:read_csv|read\.csv|read_excel|load|readRDS)\s*\(",
      "file-output":r"\b(?:ggsave|write_csv|write\.csv|saveRDS)\s*\(",
      "sampling":r"\b(?:sample_n|sample_frac|slice_sample|sample\.int)\s*\(",
    }.items() if re.search(pat,body)]
    chunks.append(dict(chunk=number,start=i+1,end=closing+1,master_phase=4<=number<=21,
      sen_reassignment=touches,flags=";".join(population_risks)))
out=root/"results/part_a/master_33831/original_rmd_audit"
out.mkdir(parents=True,exist_ok=True)
with (out/"chunk_dependency_inventory.csv").open("w",newline="") as h:
    w=csv.DictWriter(h,fieldnames=chunks[0].keys());w.writeheader();w.writerows(chunks)
print(f"[AUDIT] {len(chunks)} R chunks. Inventory: {out/'chunk_dependency_inventory.csv'}")
for c in chunks: print(f"Chunk {c['chunk']:02d}: {c['start']}-{c['end']} | master={c['master_phase']} | {c['flags']}")
if args.instrument:
    # Guardrails deliberately outside each scientific R chunk: checked only
    # after execution of chunks 4-21, and no changes inside existing chunks.
    guard='''```{r sen_master_integrity_guard_%02d, include=FALSE}
if (!exists("SEN")) stop("Master dataset SEN is missing")
if (nrow(SEN) != 33831L || dplyr::n_distinct(SEN$Run) != 33831L)
  stop("Original Part A2 master-population integrity violation after chunk %02d")
```
'''
    inserts={}
    for c in chunks:
        if 4<=c["chunk"]<=21:
            inserts[c["end"]]=guard%(c["chunk"],c["chunk"])
    result=[]
    for i,line in enumerate(lines,1):
        result.append(line)
        if i in inserts:result.append("\n"+inserts[i])
    target=root/"analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER_GUARDED.Rmd"
    target.write_text("".join(result),encoding="utf-8")
    print(f"[READY] Separate guarded copy: {target}")
