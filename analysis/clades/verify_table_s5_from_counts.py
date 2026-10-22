#!/usr/bin/env python3
"""Reanalyse published Table S5 from five-clade gene-positive counts (no AMRFinder master).
Tests mirror 15_amr_byclade.py, using ALL clades including 0% and 100% cells.
BH applies only across the 42 displayed genes; it cannot recreate an unknown
larger AMRFinder multiple-testing family or original Monte Carlo draw.
"""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, MonteCarloMethod
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[2]
CLADE={'Clade 1A':169,'Clade 1B':759,'Clade 2':894,'Clade 3':1357,'Unassigned':127}
GROUPS=list(CLADE)
SUP=dict(zip('abcdef','ᵃᵇᶜᵈᵉᶠ'))

def make_letters(groups,pmap,alpha=0.05):
    """Verified maximal-clique compact-letter display (nonsignificant pairs share a letter)."""
    def different(a,b):return pmap[tuple(sorted((a,b)))]<=alpha
    viable=[]
    for k in range(1,len(groups)+1):
        for c in itertools.combinations(groups,k):
            if all(not different(a,b) for a,b in itertools.combinations(c,2)):viable.append(c)
    maximal=[c for c in viable if not any(set(c)<set(other) for other in viable)]
    maximal.sort(key=lambda c:(-len(c),tuple(groups.index(g) for g in c)))
    if len(maximal)>26:raise ValueError('Too many CLD letters')
    assigned={g:'' for g in groups}
    for i,c in enumerate(maximal):
        for g in c:assigned[g]+=chr(97+i)
    for a,b in itertools.combinations(groups,2):
        assert bool(set(assigned[a])&set(assigned[b])) != different(a,b)
    return assigned

def count_cell(n,N,letters=''):
    pct='0%' if n==0 else '<0.1%' if 100*n/N<0.1 else f'{100*n/N:.1f}%'
    return f'{n} ({pct})'+''.join(SUP.get(k,k) for k in letters)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--counts',type=Path,default=ROOT/'data/table_s5/Table_S5_published_pure_counts.tsv')
    p.add_argument('--out',type=Path,default=ROOT/'results/table_s5_count_reanalysis')
    p.add_argument('--mc-resamples',type=int,default=100000)
    p.add_argument('--seed',type=int,default=12345)
    a=p.parse_args()
    if a.mc_resamples<9999:p.error('Monte Carlo resamples must be >=9999')
    df=pd.read_csv(a.counts,sep='\t')
    assert list(df.columns)==['Type','Gene',*GROUPS]
    assert len(df)==42 and df.Gene.is_unique and not df.isna().any().any()
    assert df.Type.value_counts().to_dict()=={'VIRULENCE':26,'AMR':13,'STRESS':3}
    assert sum(CLADE.values())==3306
    for clade,n in CLADE.items():
        x=df[clade]
        assert pd.api.types.is_integer_dtype(x) and x.between(0,n).all(),clade
    a.out.mkdir(parents=True,exist_ok=True)
    df.to_csv(a.out/'Table_S5_verified_original_counts.tsv',sep='\t',index=False)
    omnibus=[]
    for _,r in df.iterrows():
        n=int(r[GROUPS].sum())
        table=np.array([[int(r[g]),CLADE[g]-int(r[g])] for g in GROUPS])
        pval=np.nan;method='not tested';reason=''
        if n<10:reason='rarity: fewer than 10 positive genomes'
        elif n==3306:reason='invariant presence'
        else:
            result=chi2_contingency(table,correction=False)
            if (result.expected_freq<5).any():
                mc=MonteCarloMethod(n_resamples=a.mc_resamples,rng=np.random.default_rng(a.seed))
                pval=float(chi2_contingency(table,correction=False,method=mc).pvalue)
                method=f'chi-square Monte Carlo n={a.mc_resamples} seed={a.seed}'
            else:
                pval=float(result.pvalue);method='Pearson chi-square'
        omnibus.append({'Type':r.Type,'Gene':r.Gene,'Total_present':n,'Omnibus_p':pval,'FDR_BH_displayed_42_loci':np.nan,'Test_method':method,'Reason_not_tested':reason})
    om=pd.DataFrame(omnibus)
    tested=om.Omnibus_p.notna()
    if tested.any():om.loc[tested,'FDR_BH_displayed_42_loci']=multipletests(om.loc[tested,'Omnibus_p'],method='fdr_bh')[1]
    om.to_csv(a.out/'Table_S5_omnibus_BH.tsv',sep='\t',index=False)
    lookup=om.set_index('Gene');pairwise=[];summary=[]
    for _,r in df.iterrows():
        entry=lookup.loc[r.Gene];letters={g:'' for g in GROUPS}
        if pd.notna(entry.FDR_BH_displayed_42_loci) and entry.FDR_BH_displayed_42_loci<=0.05:
            prs=[]
            for g,h in itertools.combinations(GROUPS,2):
                x,y=int(r[g]),int(r[h])
                odds,raw=fisher_exact([[x,CLADE[g]-x],[y,CLADE[h]-y]],alternative='two-sided')
                prs.append({'Gene':r.Gene,'Group1':g,'Group2':h,'G1_present':x,'G1_absent':CLADE[g]-x,'G2_present':y,'G2_absent':CLADE[h]-y,'Odds_ratio':odds,'p_raw':raw})
            adjusted=multipletests([v['p_raw'] for v in prs],method='holm')[1]
            pmap={}
            for v,q in zip(prs,adjusted):
                v['p_holm']=float(q)
                pmap[tuple(sorted((v['Group1'],v['Group2'])))]=float(q)
            pairwise.extend(prs)
            letters=make_letters(GROUPS,pmap)
        row={'Type':r.Type,'Gene':r.Gene,'Total (n=3306)':int(r[GROUPS].sum())}
        for g in GROUPS:row[f'{g} (n={CLADE[g]})']=count_cell(int(r[g]),CLADE[g],letters[g])
        row.update({'Omnibus_p':entry.Omnibus_p,'FDR_BH_displayed_42_loci':entry.FDR_BH_displayed_42_loci,'Test_method':entry.Test_method})
        summary.append(row)
    pd.DataFrame(pairwise).to_csv(a.out/'Table_S5_pairwise_Fisher_Holm.tsv',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(a.out/'Table_S5_corrected_for_review.tsv',sep='\t',index=False)
    watch=['sseK3','gogB','ybtP','gyrA_S83F','mdsB','golT','sodC1','vasL','gtgA','pagO','allD']
    audit=om[om.Gene.isin(watch)]
    audit.to_csv(a.out/'Table_S5_flagged_gene_audit.tsv',sep='\t',index=False)
    with (a.out/'analysis_manifest.json').open('w') as f:json.dump({
        'source':str(a.counts),'gene_rows':len(df),'total_genomes':3306,
        'clade_sizes':CLADE,'FDR_scope':'Only the 42 published Table S5 loci',
        'monte_carlo_resamples':a.mc_resamples,'seed':a.seed,
        'comparison':'Fisher exact all 10 pairs, Holm within gene, maximal-clique CLD'
    },f,indent=2)
    print('[PASS] Published counts and clade denominators verified')
    print('[PASS] Written:',a.out)
    print(audit[['Gene','Omnibus_p','FDR_BH_displayed_42_loci']].to_string(index=False))
if __name__=='__main__':main()
