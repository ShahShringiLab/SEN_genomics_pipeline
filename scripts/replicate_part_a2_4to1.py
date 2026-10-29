#!/usr/bin/env python3
"""Reproduce a cohort close to the archived 2018-2025 QR:QS ~4:1 outcome.

This is a NEW replication scenario, not a rewrite of the historical Rmd or
the authoritative frozen selected/QC cohort. Preserve the observed source-period
sample sizes; target 80% QR among post-2017 samples IN AGGREGATE, reallocating
unattainable source-specific quotas while auditing each stratum.
"""
import csv
import os
import random
from collections import Counter
from pathlib import Path

ROOT=Path(os.environ.get("SEN_ROOT",Path(__file__).resolve().parents[1]))
MASTER=ROOT/"results/part_a/reconciled_metadata/SEN_33831_full_enriched_metadata.csv"
FROZEN=ROOT/"data/part_a/selection/Finalupdated_SEN_Genomes1.csv"
OUT=ROOT/"results/part_a/replication_4to1"
OUT.mkdir(parents=True,exist_ok=True)
SEED=123
random.seed(SEED)

def read(path):
    if not path.is_file(): raise SystemExit(f"[MISSING] {path}")
    with path.open(newline="",encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))

def period(x):
    try: y=int(float(x))
    except (ValueError,TypeError): return None
    if y<2018: return "Pre-2018"
    if y<=2020: return "2018-2020"
    if y<=2025: return "2021-2025"
    return None

full=read(MASTER)
frozen=read(FROZEN)
assert len(full)==33831 and len(frozen)==3434
assert len({r["Run"] for r in full})==33831
assert len({r["Run"] for r in frozen})==3434
for k in ("Run","Source","Collection_year","QR_genotype"):
    assert all(k in r for r in full[:1]+frozen[:1]),k

# Preserve the actual archival pre-2018 cohort, not a new historical resample.
base_ids={r["Run"] for r in frozen if period(r["Collection_year"])=="Pre-2018"}
baseline=[r for r in full if r["Run"] in base_ids]
assert len(baseline)==len(base_ids)
# Keep the achieved archival sample sizes for EACH year block/source.
targets=Counter((period(r["Collection_year"]),r["Source"]) for r in frozen
                if period(r["Collection_year"]) in ("2018-2020","2021-2025"))
assert sum(targets.values())==1965
pool={}
for r in full:
    p=period(r["Collection_year"])
    k=(p,r["Source"])
    if k in targets:
        assert str(r["QR_genotype"]).strip() in ("0","1"),(r["Run"],r["QR_genotype"])
        pool.setdefault((k,str(r["QR_genotype"]).strip()),[]).append(r)

for k,n in targets.items():
    assert len(pool.get((k,"0"),[]))+len(pool.get((k,"1"),[]))>=n,k
# Apportion QR with source-period targets; cap at physically available genomes.
qr_target=round(sum(targets.values())*0.8)
qr={}
for k,n in targets.items():
    available_qr=len(pool.get((k,"1"),[]))
    available_qs=len(pool.get((k,"0"),[]))
    lower=max(0,n-available_qs)
    upper=min(n,available_qr)
    assert lower<=upper,(k,n,available_qr,available_qs)
    qr[k]=max(lower,min(upper,round(n*0.8)))
# Reallocate deficits across feasible strata; match ~4:1 overall.
while sum(qr.values())<qr_target:
    eligible=[k for k in targets if qr[k]<min(targets[k],len(pool.get((k,"1"),[])))]
    if not eligible: raise SystemExit("Cannot meet pooled QR target with current availability")
    k=max(eligible,key=lambda z:(min(targets[z],len(pool[(z,"1")]))-qr[z],str(z)))
    qr[k]+=1
while sum(qr.values())>qr_target:
    eligible=[k for k in targets if qr[k]>max(0,targets[k]-len(pool.get((k,"0"),[])))]
    if not eligible: raise SystemExit("Cannot meet pooled QS target")
    k=max(eligible,key=lambda z:(qr[z]-max(0,targets[z]-len(pool.get((z,"0"),[]))),str(z)))
    qr[k]-=1

selected=list(baseline)
audit=[]
for k in sorted(targets):
    qn=qr[k]; sn=targets[k]-qn
    chosen=random.sample(pool.get((k,"1"),[]),qn)+random.sample(pool.get((k,"0"),[]),sn)
    selected.extend(chosen)
    audit.append({"Period":k[0],"Source":k[1],"Target_N":targets[k],
                  "Available_QR":len(pool.get((k,"1"),[])),
                  "Available_QS":len(pool.get((k,"0"),[])),
                  "Selected_QR":qn,"Selected_QS":sn,
                  "QR_pct":round(100*qn/targets[k],2)})
assert len({r["Run"] for r in selected})==len(selected)
assert len(selected)==3434,(len(selected),len(baseline))
assert sum(a["Selected_QR"] for a in audit)==qr_target
assert sum(a["Selected_QS"] for a in audit)==sum(targets.values())-qr_target
assert {r["Run"] for r in selected if r["Run"] in base_ids}==base_ids
for label,items in (("SEN_3434_replication_4to1.csv",selected),
                    ("SEN_2018_2025_4to1_source_period_audit.csv",audit)):
    with (OUT/label).open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(items[0]))
        writer.writeheader()
        writer.writerows(items)
qr_n=sum(a["Selected_QR"] for a in audit)
qs_n=sum(a["Selected_QS"] for a in audit)
print(f"[PASS] Replicated 3,434 genome candidate: {len(baseline)} pre-2018 retained; {qr_n} QR + {qs_n} QS (2018-2025); ratio {qr_n/qs_n:.3f}:1")
print(f"[OUTPUT] {OUT}")
print("[NOTE] This is a new independently sampled replication candidate, not the archived QC-adjudicated 3,434 or the original selection algorithm.")
