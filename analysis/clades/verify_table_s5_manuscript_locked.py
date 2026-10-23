#!/usr/bin/env python3
"""Reproduce the ORIGINAL manuscript Table S5 from its published count data.

This is intentionally NOT the current SEN main MonteCarloMethod(5000) analysis.

Historical provenance: the older SENBio/15_amr_byclade.py tried an unsupported
chi2_contingency(method="montecarlo",num_resamples=5000) call, caught TypeError
and fell back to Pearson chi-square without Yates correction. The user-supplied
original manuscript p-values correspond exactly to that fallback.

Keep all other original SEN analysis logic (zero/full clade exclusion, >=10
positives, BH of tested genes, pairwise Fisher/Holm, original CLD) by importing
the existing count adapter and replacing ONLY its chi-square test function.
This module never overwrites either existing results folder.
"""
import argparse
import importlib.util
import math
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import chi2_contingency

ROOT = Path(__file__).resolve().parents[2]
ADAPTER = Path(__file__).with_name("verify_table_s5_from_counts.py")
OUTPUT = ROOT / "results/table_s5_original_manuscript_locked"

# Frozen reference values transcribed from the user's original manuscript
# Table S5 on 2026-10-22. An empty field means original 'not tested'.
# These regression checks make future accidental statistical changes visible.
REFERENCE = {
    "mdsB": (0.03603943793484031, 0.04068968799094874),
    "gyrA_D87Y": (0.0, 0.0),
    "gyrA_D87N": (0.0, 0.0),
    "gyrA_D87G": (4.273602928715716e-68, 2.4929350417508346e-67),
    "gyrA_S83Y": (5.778024561936992e-06, 7.778109987222874e-06),
    "tet(A)": (1.5547663351353494e-35, 4.1859093638259405e-35),
    "qnrB19": (0.005104235788510471, 0.0059549417532622165),
    "blaTEM-1": (2.670518186953949e-41, 8.497103322126202e-41),
    "aph(3'')-Ib": (1.4669122646649866e-70, 1.0268385852654906e-69),
    "sul2": (5.536897668799093e-62, 2.4223927300996034e-61),
    "aph(6)-Id": (6.331196935850621e-52, 2.2159189275477177e-51),
    "blaTEM": (6.045634326669722e-19, 1.3224825089590018e-18),
    "gyrA_S83F": (0.2798030414215197, 0.2967608015076724),
    "golT": (None, None),
    "merR": (3.216988208123137e-11, 5.361647013538562e-11),
    "qacEdelta1": (3.57437416235274e-07, 5.212628986764412e-07),
    "rhuM": (0.004368991432013139, 0.00527292069380896),
    "sefA": (0.00031638622943526407, 0.0003954827867940801),
    "allD": (None, None),
    "avrA": (0.0001454175800308172, 0.00018850427041031858),
    "sinH": (0.44562619987526764, 0.44562619987526764),
    "pagO": (0.11984957262591724, 0.13108547005959698),
    "sodC1": (None, None),
    "vasL": (None, None),
    "gtgA": (None, None),
    "sseI": (5.092591859046613e-08, 7.749596307244845e-08),
    "bfpA": (5.509917090440559e-245, 6.428236605513986e-244),
    "sspH2": (3.341928135452421e-41, 9.74729039506956e-41),
    "spvB": (5.791489906910648e-13, 1.0668534039045932e-12),
    "mltE": (2.416911738507559e-12, 4.229595542388229e-12),
    "rcK": (3.394870767639922e-78, 2.970511921684932e-77),
    "gtrC1": (1.656612290483799e-08, 2.6355195530424075e-08),
    "rfc": (1.6881617378151657e-06, 2.3634264329412318e-06),
    "sty2629": (1.0701759573883427e-16, 2.203303441681882e-16),
    "spvD": (2.167233825807325e-13, 4.2140657724031323e-13),
    "sseK3": (0.32601598275860444, 0.3356046881338575),
    "int": (1.8556658364775664e-63, 9.278329182387832e-63),
    "xis": (1.4137544648621037e-60, 5.497934030019293e-60),
    "safB": (8.80488393860053e-24, 2.0544729190067906e-23),
    "traT": (7.265462251696637e-31, 1.816365562924159e-30),
    "gogB": (None, None),
    "ybtP": (None, None),
}


def original_legacy_fallback(table):
    """The *actual* result of the original legacy chi2 Monte Carlo fallback."""
    result = chi2_contingency(table, correction=False)
    return float(result.pvalue), "chi-square", float(result.statistic), int(result.dof)


def assert_original_reference(outdir: Path) -> None:
    results = outdir / "Split_Clade1A1B/With_Unassigned/AMRFinder_Omnibus_Stats.tsv"
    reported = pd.read_csv(results, sep="\t").set_index("Element symbol")
    if set(reported.index) != set(REFERENCE):
        raise AssertionError(f"Reference gene set mismatch: {set(REFERENCE)^set(reported.index)}")
    mismatches = []
    for gene, (exp_p, exp_fdr) in REFERENCE.items():
        obs = reported.loc[gene]
        for label, actual, expected in (
            ("p-value", obs["Stats p-value"], exp_p),
            ("FDR", obs["Stats FDR"], exp_fdr),
        ):
            if expected is None:
                ok = pd.isna(actual)
            else:
                ok = pd.notna(actual) and math.isclose(
                    float(actual), expected, rel_tol=2e-11, abs_tol=0
                )
            if not ok:
                mismatches.append((gene, label, expected, actual))
    # In addition, verify that the original significance alphabets survive.
    manuscript = pd.read_excel(
        outdir/"Split_Clade1A1B/With_Unassigned/Table_S5_Manuscript_Ready.xlsx",
        sheet_name="Manuscript_Table_S5",
    ).set_index("Gene")
    checks = {
        ("mdsB","Clade 3 (n=1357)"): "ᵃ",
        ("mdsB","Unassigned (n=127)"): "ᵃ",
        ("gyrA_D87Y","Clade 1A (n=169)"): "ᵃ",
        ("gyrA_D87Y","Clade 3 (n=1357)"): "ᵈ",
        ("sseK3","Clade 1A (n=169)"): None,
    }
    for (gene, col), expected in checks.items():
        value = str(manuscript.loc[gene,col])
        if expected is None:
            ok = not any(c in value for c in "ᵃᵇᶜᵈᵉᶠ")
        else:
            ok = expected in value
        if not ok:
            mismatches.append((gene,col,expected,value))
    if mismatches:
        for entry in mismatches:
            print("[MISMATCH]",entry,file=sys.stderr)
        raise AssertionError(f"{len(mismatches)} historical manuscript discrepancies")
    print(f"[PASS] ALL {len(REFERENCE)} genes reproduce the original manuscript p-values and FDR.")
    print("[PASS] Checked manuscript superscript letters against original.")
    print("[OUTPUT]",outdir)


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts",type=Path,
        default=ROOT/"data/table_s5/Table_S5_published_pure_counts.tsv")
    parser.add_argument("--out",type=Path,default=OUTPUT)
    args = parser.parse_args()
    if not ADAPTER.is_file():
        raise SystemExit(f"Missing base adapter: {ADAPTER}")
    spec = importlib.util.spec_from_file_location("sen_counts_main_adapter", ADAPTER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Retain original pairwise, BH and letter source; only restore the
    # historical Pearson chi-square fallback used to create this manuscript.
    module.chi2_with_optional_montecarlo = original_legacy_fallback
    old_argv = sys.argv
    try:
        sys.argv = [str(ADAPTER),"--counts",str(args.counts),"--out",str(args.out)]
        module.run()
    finally:
        sys.argv = old_argv
    assert_original_reference(args.out)
    (args.out/"STATISTICAL_PROVENANCE.txt").write_text(
        "ORIGINAL manuscript-compatible SEN Table S5.\n"
        "Reference: user-provided 42-row Table S5 including all original p-values/FDR.\n"
        "Original legacy Monte Carlo call failed and fell back to Pearson "
        "chi-square (correction=False). This script reproduces that behavior.\n"
        "Not equivalent to current SEN main seeded Monte Carlo implementation.\n"
        "All other statistical rules come from the existing count-only adapter.\n"
        "No sample-level AMRFinder annotations are reconstructed.\n"
    )

if __name__ == "__main__":
    run()
