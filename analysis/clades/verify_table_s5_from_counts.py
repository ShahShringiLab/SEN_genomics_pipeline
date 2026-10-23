#!/usr/bin/env python3
"""Historical ORIGINAL 15_amr_byclade.py statistics, applied to published Table S5 counts.

No corrections or statistical substitutions. Exact original function bodies
are copied from ajulojays/SENBio-Final main:15_amr_byclade.py.

Limitations of a count-only subset: original AMRFinder master contains additional
loci that change BH FDR across the complete testing family; original code's
legacy Monte Carlo API may TypeError in modern SciPy and silently use Pearson.
These caveats are printed/reported; the original behavior is NOT modified.
"""
import argparse
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[2]
clade_order=["Clade 1A","Clade 1B","Clade 2","Clade 3","Unassigned"]
clade_sizes={"Clade 1A":169,"Clade 1B":759,"Clade 2":894,"Clade 3":1357,"Unassigned":127}
MIN_TOTAL_PRESENT=10
ALPHA=0.05
SUPERSCRIPT={"a":"ᵃ","b":"ᵇ","c":"ᶜ","d":"ᵈ","e":"ᵉ","f":"ᶠ","g":"ᵍ","h":"ʰ","i":"ᶦ"}

# Verbatim statistical functions from the original pre-correction script.
def chi2_with_optional_montecarlo(table: np.ndarray):
    chi2, p, dof, expected = chi2_contingency(table, correction=False)
    method = "chi-square"
    if np.any(expected < 5):
        try:
            chi2_mc, p_mc, dof_mc, _ = chi2_contingency(
                table, correction=False, method="montecarlo", num_resamples=5000
            )
            return float(p_mc), "chi-square (montecarlo)", float(chi2_mc), int(dof_mc)
        except TypeError:
            pass
        except Exception:
            pass
    return float(p), method, float(chi2), int(dof)

def odds_ratio_with_haldene(a, b, c, d, add=0.5):
    return ((a + add) / (b + add)) / ((c + add) / (d + add))

def make_cld_letters(groups, p_adj_map, alpha=0.05):
    if not p_adj_map or len(groups) == 0:
        return {g: "" for g in groups}

    def different(g1, g2):
        if g1 == g2:
            return False
        key = tuple(sorted((g1, g2)))
        p = p_adj_map.get(key, np.nan)
        return (pd.notna(p) and p <= alpha)

    letter_pool = list("abcdefghijklmnopqrstuvwxyz")
    assign = {g: set() for g in groups}
    remaining = set(groups)

    while remaining and letter_pool:
        L = letter_pool.pop(0)
        block = []
        for g in groups:
            if g not in remaining:
                continue
            ok = True
            for h in block:
                if different(g, h):
                    ok = False
                    break
            if ok:
                block.append(g)

        if not block:
            g = next(iter(remaining))
            block = [g]

        for g in block:
            assign[g].add(L)
            if g in remaining:
                remaining.remove(g)

    return {g: "".join(sorted(assign[g])) for g in groups}


def pct_format(count,denom):
    if denom<=0:return ""
    if count==0:return "0%"
    pct=100.0*count/denom
    if 0<pct<0.1:return "<0.1%"
    if abs(pct-round(pct))<1e-12:return f"{int(round(pct))}%"
    return f"{pct:.1f}%"

def format_cell(count,denom,sup_letters=""):
    p=pct_format(count,denom)
    if p=="":return ""
    if sup_letters:return f"{count} ({p}){''.join(SUPERSCRIPT.get(ch,ch) for ch in sup_letters)}"
    return f"{count} ({p})"

def run():
    parser=argparse.ArgumentParser(description="Historical exact 15_amr_byclade.py statistical logic applied to Table S5 counts")
    parser.add_argument("--counts",type=Path,default=ROOT/"data/table_s5/Table_S5_published_pure_counts.tsv")
    parser.add_argument("--out",type=Path,default=ROOT/"results/table_s5_historical_replication")
    args=parser.parse_args()
    df=pd.read_csv(args.counts,sep="\t")
    assert list(df.columns)==["Type","Gene",*clade_order]
    assert len(df)==42 and df["Gene"].is_unique
    for c in clade_order:
        assert pd.api.types.is_integer_dtype(df[c]) and df[c].between(0,clade_sizes[c]).all()
    args.out.mkdir(parents=True,exist_ok=True)
    df.to_csv(args.out/"Table_S5_original_pure_counts.tsv",sep="\t",index=False)
    records=[]
    used_clades_by_element={}
    for _,r in df.iterrows():
        symbol=r["Gene"]
        pres={c:int(r[c]) for c in clade_order}
        ns={c:clade_sizes[c] for c in clade_order}
        absn={c:ns[c]-pres[c] for c in clade_order}
        tot_pres=sum(pres.values())
        # Original rule - DO NOT change or include 0%/100% clades.
        variable_clades=[c for c in clade_order if ns[c]>0 and (0<pres[c]<ns[c])]
        used_clades_by_element[symbol]=variable_clades
        pval=np.nan;method="";chi2=np.nan;dof=np.nan;reason=""
        if tot_pres<MIN_TOTAL_PRESENT:
            reason=f"rare (<{MIN_TOTAL_PRESENT} total present)"
        elif len(variable_clades)<2:
            reason="insufficient variable clades after zero/full exclusion"
        else:
            table=np.array([[pres[c],absn[c]] for c in variable_clades],dtype=int)
            try:
                pval,method,chi2,dof=chi2_with_optional_montecarlo(table)
            except Exception as e:
                reason=f"chi2 failed: {e}"
        records.append({"Element symbol":symbol,"Type":r["Type"],"Total_Present":tot_pres,
                        "Tested_Clades":",".join(variable_clades),"Test_Method":method,
                        "Chi2":chi2,"DoF":dof,"Stats p-value":pval,
                        "Stats FDR":np.nan,"Reason_Not_Tested":reason})
    omnibus_df=pd.DataFrame(records)
    p_num=pd.to_numeric(omnibus_df["Stats p-value"],errors="coerce")
    tested_mask=p_num.notna()
    if tested_mask.sum()>0:
        _,fdr_vals,_,_=multipletests(p_num.loc[tested_mask].values.astype(float),method="fdr_bh")
        omnibus_df.loc[tested_mask,"Stats FDR"]=fdr_vals
    omnibus_df.to_csv(args.out/"AMRFinder_Omnibus_Stats.tsv",sep="\t",index=False)
    fdr_num=pd.to_numeric(omnibus_df["Stats FDR"],errors="coerce")
    sig_elements=omnibus_df.loc[tested_mask&(fdr_num<=ALPHA),"Element symbol"].tolist()
    pairwise_rows=[]
    cld_letters_by_element={e:{c:"" for c in clade_order} for e in df.Gene}
    for symbol in sig_elements:
        row=df.loc[df.Gene==symbol].iloc[0]
        pres={c:int(row[c]) for c in clade_order}
        ns={c:clade_sizes[c] for c in clade_order}
        absn={c:ns[c]-pres[c] for c in clade_order}
        variable_clades=used_clades_by_element.get(symbol,[])
        if len(variable_clades)<2:continue
        pairs=list(itertools.combinations(variable_clades,2))
        pvals=[];tmp_rows=[]
        for a,b in pairs:
            a1,a0=pres[a],absn[a];b1,b0=pres[b],absn[b]
            _,p_raw=fisher_exact([[a1,a0],[b1,b0]],alternative="two-sided")
            or05=odds_ratio_with_haldene(a1,a0,b1,b0,add=0.5)
            tmp_rows.append({"Element symbol":symbol,"Group1":a,"Group2":b,
                             "G1_present":a1,"G1_absent":a0,"G2_present":b1,
                             "G2_absent":b0,"OR_0.5":or05,"p_raw":float(p_raw),"p_holm":np.nan})
            pvals.append(float(p_raw))
        if pvals:
            _,p_holm,_,_=multipletests(pvals,method="holm")
            for i in range(len(tmp_rows)):tmp_rows[i]["p_holm"]=float(p_holm[i])
        pairwise_rows.extend(tmp_rows)
        p_adj_map={}
        for rr in tmp_rows:
            key=tuple(sorted((rr["Group1"],rr["Group2"])))
            p_adj_map[key]=rr["p_holm"]
        cld_var=make_cld_letters(variable_clades,p_adj_map,alpha=ALPHA)
        for c in variable_clades:cld_letters_by_element[symbol][c]=cld_var.get(c,"")
    pd.DataFrame(pairwise_rows).to_csv(args.out/"AMRFinder_Pairwise_Fisher.tsv",sep="\t",index=False)
    p_map=dict(zip(omnibus_df["Element symbol"],omnibus_df["Stats p-value"]))
    f_map=dict(zip(omnibus_df["Element symbol"],omnibus_df["Stats FDR"]))
    summary=[]
    for _,r in df.iterrows():
        symbol=r.Gene
        item={"Type":r.Type,"Element symbol":symbol,"Total (n)":int(r[clade_order].sum())}
        letters=cld_letters_by_element[symbol]
        for c in clade_order:
            item[f"{c} (n={clade_sizes[c]})"]=format_cell(int(r[c]),clade_sizes[c],letters[c])
        item["Stats p-value"]="" if pd.isna(p_map[symbol]) else float(p_map[symbol])
        item["Stats FDR"]="" if pd.isna(f_map[symbol]) else float(f_map[symbol])
        summary.append(item)
    summary_df=pd.DataFrame(summary)
    summary_df.to_csv(args.out/"AMRFinder_ByClade_Summary.tsv",sep="\t",index=False)
    print("[PASS] Exact historical statistical selection, testing, Holm and greedy-letter logic run on 42 published genes.")
    print("[CAVEAT] Historical Monte Carlo call uses an unsupported legacy signature on SciPy 1.15.3; on TypeError it falls back to asymptotic chi-square, as in the original script.")
    print("[CAVEAT] BH FDR covers only the 42 displayed genes, not all original AMRFinder elements.")
    print("[OUTPUT]",args.out)
    print(omnibus_df.loc[omnibus_df["Element symbol"].isin(["sseK3","gogB","ybtP","gyrA_S83F"])].to_string(index=False))
if __name__=="__main__":run()
