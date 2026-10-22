#!/usr/bin/env python3
"""Lightweight regression tests for AMRFinder-by-clade Table S5 logic.

Extracts the two pure statistic functions from the production script via AST,
so importing the analysis does not try to open the 133k-row master file.
"""
import ast
import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, MonteCarloMethod

root=Path(__file__).resolve().parents[1]
script=(root/"analysis/clades/15_amr_byclade.py").read_text()
parsed=ast.parse(script)
functions=[n for n in parsed.body if isinstance(n,ast.FunctionDef) and
           n.name in {"make_cld_letters","chi2_with_optional_montecarlo"}]
assert len(functions)==2
env={"np":np,"pd":pd,"itertools":itertools,"chi2_contingency":chi2_contingency,
     "MonteCarloMethod":MonteCarloMethod,
     "args":argparse.Namespace(mc_resamples=10000)}
exec(compile(ast.Module(body=functions,type_ignores=[]),filename=str(root),mode="exec"),env)
cld=env["make_cld_letters"]
chi=env["chi2_with_optional_montecarlo"]

groups=["Clade 1A","Clade 1B","Clade 2","Clade 3","Unassigned"]
def pair_map(val):
    return {tuple(sorted((a,b))):val for a,b in itertools.combinations(groups,2)}
all_same=cld(groups,pair_map(1.0))
assert len(set(all_same.values()))==1, all_same
all_different=cld(groups,pair_map(1e-10))
assert len(set(all_different.values()))==5, all_different
assert all(not(set(all_different[a])&set(all_different[b]))
           for a,b in itertools.combinations(groups,2))

p=pair_map(0.5)
p[tuple(sorted(("Clade 1A","Clade 2")))]=0.001
letters=cld(groups,p)
assert set(letters["Clade 1A"]) & set(letters["Clade 1B"])
assert set(letters["Clade 1B"]) & set(letters["Clade 2"])
assert not(set(letters["Clade 1A"]) & set(letters["Clade 2"])),letters

sizes=np.array([169,759,894,1357,127],dtype=int)
for symbol,present in {
    "sseK3":[10,70,0,0,13],
    "gogB":[0,0,0,0,22],
    "ybtP":[0,0,0,0,10],
    "gyrA_S83F":[0,3,6,0,2],
}.items():
    present=np.array(present,dtype=int)
    table=np.column_stack((present,sizes-present))
    pval,method,stat,dof=chi(table)
    assert np.isfinite(pval) and pval<0.05,(symbol,pval,method)
    assert table.shape==(5,2),symbol
    print(f"[PASS] {symbol}: 5-clade table, n={sum(present)}, p={pval:.7g}, method={method}")
assert "(0 < pres[c] < ns[c])" not in script
assert "Tested_Clades" in script and "tested_clades" in script
assert "if args.stats_only:" in script
assert "must have one explicit Sample column" in script
print("[PASS] All pairwise CLD and fixed-zero/full five-clade omnibus regression tests passed")
