#!/usr/bin/env python3
import os
import re
from pathlib import Path
import itertools
import warnings
import argparse
import json
from datetime import datetime, timezone

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import chi2_contingency, fisher_exact, MonteCarloMethod
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
# CONFIG
# =============================================================================
REPO_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = str(Path(os.environ.get("SEN_ROOT", REPO_ROOT)))

parser = argparse.ArgumentParser(description="AMRFinder clade stratification and audited Table S5 statistics")
parser.add_argument("--amr-master", default=os.environ.get("SEN_AMRFINDER_MASTER", os.path.join(BASE_DIR, "AMRFinderplus", "master_AMRFinder_report.tsv")))
parser.add_argument("--clade-metadata", default=os.environ.get("SEN_CLADE_METADATA", os.path.join(BASE_DIR, "metadata", "final_clade_metadata.tsv")))
parser.add_argument("--scheme", choices=["Combined_Clade1", "Split_Clade1A1B"], help="Only run a single scheme; default both")
parser.add_argument("--mode", choices=["With_Unassigned", "Without_Unassigned"], help="Only run a single mode; default both")
parser.add_argument("--stats-only", action="store_true", help="Skip extremely large heatmaps but keep TSV and Excel statistics")
parser.add_argument("--mc-resamples", type=int, default=100000, help="Monte Carlo resamples for sparse omnibus tests (default 100000)")
args = parser.parse_args()
if args.mc_resamples < 9999:
    parser.error("--mc-resamples must be >= 9999 for scientifically meaningful sparse-test resolution")
amr_report_path = args.amr_master

metadata_path_candidates = [
    args.clade_metadata,
    os.path.join(BASE_DIR, "iqtree_final", "ITOL", "Clade_metadata.txt"),
    os.path.join(BASE_DIR, "ITOL", "Clade_metadata.txt"),
    os.path.join(BASE_DIR, "iqtree", "ITOL", "Clade_metadata.txt"),
]
metadata_path = next((p for p in metadata_path_candidates if os.path.exists(p)), None)
if metadata_path is None:
    raise SystemExit("[ERROR] Could not find Clade_metadata.txt in expected locations.")

ROOT_OUT_DIR = os.path.join(BASE_DIR, "iqtree_final", "AMRFinder_clade")
os.makedirs(ROOT_OUT_DIR, exist_ok=True)

# Presence rule
MIN_COVERAGE = 80.0
MIN_IDENTITY = 90.0
MIN_TOTAL_PRESENT = 10
ALPHA = 0.05

# Remove SMOKE rows entirely
EXCLUDE_SMOKE_PREFIX = "_SMOKE_"

# Base clade orders
BASE_CLADE_ORDER_SPLIT = ["Clade 1A", "Clade 1B", "Clade 2", "Clade 3", "Unassigned"]
BASE_CLADE_ORDER_COMBINED = ["Clade 1", "Clade 2", "Clade 3", "Unassigned"]

# Clade colors
CLADE_PALETTE_SPLIT = {
    "Clade 1A": "#FF7F0E",
    "Clade 1B": "#FF0000",
    "Clade 2": "#0000FF",
    "Clade 3": "#008000",
    "Unassigned": "#808080",
}
CLADE_PALETTE_COMBINED = {
    "Clade 1": "#FF0000",
    "Clade 2": "#0000FF",
    "Clade 3": "#008000",
    "Unassigned": "#808080",
}

# Presence/absence heatmap colors
my_cmap = ListedColormap(["#FFFFFF", "#FF0000"])

# =============================================================================
# HELPERS
# =============================================================================
SUPERSCRIPT = {
    "a": "ᵃ", "b": "ᵇ", "c": "ᶜ", "d": "ᵈ", "e": "ᵉ", "f": "ᶠ", "g": "ᵍ", "h": "ʰ", "i": "ᶦ",
    "j": "ʲ", "k": "ᵏ", "l": "ˡ", "m": "ᵐ", "n": "ⁿ", "o": "ᵒ", "p": "ᵖ", "r": "ʳ", "s": "ˢ",
    "t": "ᵗ", "u": "ᵘ", "v": "ᵛ", "w": "ʷ", "x": "ˣ", "y": "ʸ", "z": "ᶻ"
}

def letters_to_sup(letters: str) -> str:
    return "".join(SUPERSCRIPT.get(ch, ch) for ch in letters)

def normalize_sample_id(x: str) -> str:
    if pd.isna(x):
        return x
    s = str(x).strip()
    s = os.path.basename(s)
    s = re.sub(r"\.tsv$", "", s, flags=re.IGNORECASE)
    s = s.replace("_trimmed", "").replace("_pure", "")
    s = s.strip("_")
    m = re.search(r"(SRR\d+)", s)
    return m.group(1) if m else s

def safe_type_palette(type_values):
    types = sorted(list(pd.unique(pd.Series(type_values).dropna())))
    pal = sns.color_palette("tab20", n_colors=max(1, len(types)))
    return {t: pal[i % len(pal)] for i, t in enumerate(types)}

def cluster_columns_within_group(binary_df: pd.DataFrame) -> list:
    cols = list(binary_df.columns)
    if len(cols) <= 2:
        return cols
    X = binary_df[cols].T.values.astype(float)
    try:
        Z = linkage(X, method="average", metric="euclidean")
        order = leaves_list(Z)
        return [cols[i] for i in order]
    except Exception:
        return cols

def pct_format(count: int, denom: int) -> str:
    if denom <= 0:
        return ""
    if count == 0:
        return "0%"
    pct = 100.0 * count / denom
    if 0 < pct < 0.1:
        return "<0.1%"
    if abs(pct - round(pct)) < 1e-12:
        return f"{int(round(pct))}%"
    return f"{pct:.1f}%"

def format_cell(count: int, denom: int, sup_letters: str = "") -> str:
    p = pct_format(count, denom)
    if p == "":
        return ""
    if sup_letters:
        return f"{count} ({p}){letters_to_sup(sup_letters)}"
    return f"{count} ({p})"

def tsv_integrity_check(path: str):
    bad = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
        expected_nf = header.count("\t") + 1
        for i, line in enumerate(fh, start=2):
            nf = line.count("\t") + 1
            if nf != expected_nf:
                bad.append({"line": i, "fields": nf, "expected": expected_nf})
    return expected_nf, bad

def read_master_tsv_robust(path: str, badlines_out: str):
    expected_nf, bad = tsv_integrity_check(path)
    print(f"[INFO] master header fields = {expected_nf}")

    if bad:
        print(f"[FAIL] Found {len(bad)} malformed lines in master TSV.")
        pd.DataFrame(bad).to_csv(badlines_out, sep="\t", index=False)
        raise SystemExit(f"[FAIL] Refusing to silently omit AMRFinder rows. Inspect {badlines_out}")

    return pd.read_csv(path, sep="\t", low_memory=False)

def chi2_with_optional_montecarlo(table: np.ndarray):
    result = chi2_contingency(table, correction=False)
    expected = result.expected_freq

    if np.any(expected < 5):
        # Fixed seed + fixed number of resamples makes sparse-count handling
        # reproducible. SciPy >=1.15 requires a MonteCarloMethod instance.
        mc = MonteCarloMethod(
            n_resamples=args.mc_resamples,
            rng=np.random.default_rng(12345),
        )
        mc_result = chi2_contingency(
            table,
            correction=False,
            method=mc,
        )
        return (
            float(mc_result.pvalue),
            f"chi-square (montecarlo, {args.mc_resamples} resamples, seed=12345)",
            float(mc_result.statistic),
            np.nan,
        )

    return (
        float(result.pvalue),
        "chi-square",
        float(result.statistic),
        int(result.dof),
    )

def odds_ratio_with_haldene(a, b, c, d, add=0.5):
    return ((a + add) / (b + add)) / ((c + add) / (d + add))

def make_cld_letters(groups, p_adj_map, alpha=0.05):
    """Compact letter display: every non-significant pair shares a letter,
    every significant pair shares none. Handles nontransitive Holm results
    (e.g. A~B, B~C, yet A differs from C) via overlapping maximal cliques.
    """
    groups = list(groups)
    if not groups:
        return {}
    def significant(a,b):
        p = p_adj_map.get(tuple(sorted((a,b))), np.nan)
        if not np.isfinite(p):
            raise ValueError(f"Missing Holm-adjusted comparison: {a} vs {b}")
        return bool(p <= alpha)
    nonsig_cliques = []
    for size in range(1,len(groups)+1):
        for clique in itertools.combinations(groups,size):
            if all(not significant(a,b) for a,b in itertools.combinations(clique,2)):
                nonsig_cliques.append(clique)
    maximal = [c for c in nonsig_cliques
               if not any(set(c)<set(other) for other in nonsig_cliques)]
    maximal.sort(key=lambda c:(-len(c),tuple(groups.index(x) for x in c)))
    pool = "abcdefghijklmnopqrstuvwxyz"
    if len(maximal)>len(pool):
        raise ValueError("Too many CLD maximal cliques for the letter pool")
    assigned = {g:"" for g in groups}
    for idx,clique in enumerate(maximal):
        for g in clique:
            assigned[g] += pool[idx]
    for a,b in itertools.combinations(groups,2):
        shares=bool(set(assigned[a]) & set(assigned[b]))
        assert shares != significant(a,b), f"CLD inconsistent: {a}/{b}"
    return assigned

def recode_lineages(meta_df: pd.DataFrame, scheme: str) -> pd.DataFrame:
    out = meta_df.copy()
    if scheme == "Combined_Clade1":
        out["Lineage"] = out["Lineage"].replace({
            "Clade 1A": "Clade 1",
            "Clade 1B": "Clade 1",
        })
    elif scheme == "Split_Clade1A1B":
        pass
    else:
        raise ValueError(f"Unknown scheme: {scheme}")
    return out

def get_scheme_config(scheme: str, mode: str):
    include_unassigned = (mode == "With_Unassigned")

    if scheme == "Combined_Clade1":
        clade_order = ["Clade 1", "Clade 2", "Clade 3"]
        palette = CLADE_PALETTE_COMBINED.copy()
    elif scheme == "Split_Clade1A1B":
        clade_order = ["Clade 1A", "Clade 1B", "Clade 2", "Clade 3"]
        palette = CLADE_PALETTE_SPLIT.copy()
    else:
        raise ValueError(f"Unknown scheme: {scheme}")

    if include_unassigned:
        clade_order = clade_order + ["Unassigned"]

    return clade_order, palette, include_unassigned

def build_output_paths(base_out: str):
    return {
        "excel": os.path.join(base_out, "AMRFinder_ByClade_Summary.xlsx"),
        "table": os.path.join(base_out, "AMRFinder_ByClade_ElementCounts.tsv"),
        "matrix": os.path.join(base_out, "AMRFinder_PresenceAbsence_Matrix.tsv"),
        "badlines": os.path.join(base_out, "master_bad_lines.tsv"),
        "heat_png": os.path.join(base_out, "AMRFinder_Heatmap_ByType.png"),
        "heat_pdf": os.path.join(base_out, "AMRFinder_Heatmap_ByType.pdf"),
        "omnibus": os.path.join(base_out, "AMRFinder_Omnibus_Stats.tsv"),
        "pairwise": os.path.join(base_out, "AMRFinder_Pairwise_Fisher.tsv"),
        "summary": os.path.join(base_out, "AMRFinder_ByClade_Summary.tsv"),
    }

# =============================================================================
# LOAD MASTER DATA ONCE
# =============================================================================
print("📂 Loading data...")

base_out_for_badlines = os.path.join(ROOT_OUT_DIR, "_tmp")
os.makedirs(base_out_for_badlines, exist_ok=True)
tmp_badlines = os.path.join(base_out_for_badlines, "master_bad_lines.tsv")

if not os.path.exists(amr_report_path):
    raise SystemExit(f"[ERROR] Missing: {amr_report_path}")
if not os.path.exists(metadata_path):
    raise SystemExit(f"[ERROR] Missing: {metadata_path}")

amr_df = read_master_tsv_robust(amr_report_path, tmp_badlines)
meta_df_master = pd.read_csv(metadata_path, sep="\t", low_memory=False)

if "Sample" not in meta_df_master.columns or "Lineage" not in meta_df_master.columns:
    raise SystemExit("[ERROR] Clade_metadata.txt must have columns: Sample, Lineage")

amr_df.rename(columns={amr_df.columns[0]: "Sample"} if amr_df.columns[0] != "Sample" else {}, inplace=True)

amr_df["Sample_raw"] = amr_df["Sample"].astype(str)
amr_df = amr_df[~amr_df["Sample_raw"].str.startswith(EXCLUDE_SMOKE_PREFIX)].copy()

amr_df["Sample"] = amr_df["Sample_raw"].map(normalize_sample_id)
meta_df_master["Sample"] = meta_df_master["Sample"].map(normalize_sample_id)

meta_df_master = meta_df_master[meta_df_master["Lineage"] != "Reference"].copy()

required_lineages = {"Clade 1A", "Clade 1B", "Clade 2", "Clade 3", "Unassigned"}
unknown_lineages = set(meta_df_master["Lineage"].dropna().unique()) - required_lineages
if unknown_lineages or meta_df_master["Lineage"].isna().any():
    raise SystemExit(f"[FAIL] Unexpected clade labels (not silently Unassigned): {unknown_lineages}")
if meta_df_master["Sample"].isna().any() or meta_df_master["Sample"].duplicated().any():
    raise SystemExit("[FAIL] Duplicate or missing sample identifiers in clade assignments")
if len(meta_df_master) != 3306:
    raise SystemExit(f"[FAIL] Expected 3306 unique assigned samples, got {len(meta_df_master)}")

req_cols = {"Sample", "Type", "Element symbol", "Element name"}
missing = req_cols - set(amr_df.columns)
if missing:
    raise SystemExit(f"[ERROR] master_AMRFinder_report.tsv missing columns: {missing}")

# robust coverage/identity detection
colmap = {str(c).strip().lower(): c for c in amr_df.columns}

coverage_col = None
for cand in [
    "% coverage of reference",
    "% coverage",
    "%coverage",
    "coverage"
]:
    if cand in colmap:
        coverage_col = colmap[cand]
        break

identity_col = None
for cand in [
    "% identity to reference",
    "% identity",
    "%identity",
    "identity"
]:
    if cand in colmap:
        identity_col = colmap[cand]
        break

before = len(amr_df)

if coverage_col is not None:
    amr_df[coverage_col] = pd.to_numeric(amr_df[coverage_col], errors="coerce")
    amr_df = amr_df[amr_df[coverage_col] >= MIN_COVERAGE].copy()
else:
    print("[WARN] No coverage column found in master_AMRFinder_report.tsv; skipping coverage filter.")

if identity_col is not None:
    amr_df[identity_col] = pd.to_numeric(amr_df[identity_col], errors="coerce")
    amr_df = amr_df[amr_df[identity_col] >= MIN_IDENTITY].copy()
else:
    print("[WARN] No identity column found in master_AMRFinder_report.tsv; skipping identity filter.")

after = len(amr_df)
print(f"[INFO] Applied filters: coverage >= {MIN_COVERAGE} and identity >= {MIN_IDENTITY}: {before} -> {after}")
print(f"[INFO] Coverage column used: {coverage_col}")
print(f"[INFO] Identity column used: {identity_col}")

print(f"[INFO] Raw metadata lineages found: {sorted(meta_df_master['Lineage'].dropna().unique())}")

# =============================================================================
# RUN BOTH SCHEMES × MODES
# =============================================================================
print(f"[AUDIT] Inputs: AMRFinder={amr_report_path}; clades={metadata_path}")
print("[AUDIT] Testing uses ALL nonempty clades (including 0% and 100% cells).")
print(f"[AUDIT] MC resamples={args.mc_resamples}, seed=12345; BH across all tested elements; pairwise Holm.")
for scheme in ([args.scheme] if args.scheme else ["Combined_Clade1", "Split_Clade1A1B"]):
    for mode in ([args.mode] if args.mode else ["With_Unassigned", "Without_Unassigned"]):
        print("\n" + "=" * 80)
        print(f"🚀 Processing scheme={scheme} | mode={mode}")
        print("=" * 80)

        output_dir = os.path.join(ROOT_OUT_DIR, scheme, mode)
        os.makedirs(output_dir, exist_ok=True)
        out = build_output_paths(output_dir)

        meta_df = recode_lineages(meta_df_master, scheme)
        clade_order, clade_palette, include_unassigned = get_scheme_config(scheme, mode)

        if not include_unassigned:
            meta_df = meta_df[meta_df["Lineage"] != "Unassigned"].copy()

        meta_df["Lineage"] = meta_df["Lineage"].where(meta_df["Lineage"].isin(clade_order), other="Unassigned")
        meta_df["Lineage"] = pd.Categorical(meta_df["Lineage"], categories=clade_order, ordered=True)
        meta_df = meta_df.sort_values(["Lineage", "Sample"]).reset_index(drop=True)

        allowed_samples = meta_df["Sample"].tolist()
        print(f"[INFO] Metadata samples (non-Reference, mode filtered): {len(allowed_samples)}")
        print(f"[INFO] Scheme lineages in use: {clade_order}")

        amr_curr = amr_df[amr_df["Sample"].isin(set(allowed_samples))].copy()
        print(f"[INFO] AMR rows after filters={len(amr_curr)}")

        catalog = (
            amr_curr.groupby(["Element symbol", "Type", "Element name"], dropna=False)
                  .size()
                  .reset_index(name="n")
                  .sort_values(["Element symbol", "n"], ascending=[True, False])
                  .drop_duplicates(subset=["Element symbol"], keep="first")
                  .drop(columns=["n"])
                  .set_index("Element symbol")
        )
        sym_to_type = catalog["Type"].to_dict()
        sym_to_name = catalog["Element name"].to_dict()

        print("🧱 Building presence/absence matrix...")
        pa = amr_curr.pivot_table(
            index="Sample",
            columns="Element symbol",
            aggfunc="size",
            fill_value=0
        )
        pa = (pa > 0).astype(int)

        missing_samples = meta_df.loc[~meta_df["Sample"].isin(pa.index), "Sample"].tolist()
        if missing_samples:
            miss_path = os.path.join(output_dir, "samples_in_metadata_missing_in_master.tsv")
            pd.DataFrame({"Sample": missing_samples}).to_csv(miss_path, sep="\t", index=False)
            print(f"[WARN] {len(missing_samples)} metadata samples missing in master hits; added as all-zero rows. Saved: {miss_path}")

        pa = pa.reindex(meta_df["Sample"], fill_value=0)
        pa.to_csv(out["matrix"], sep="\t")
        print(f"[OK] Matrix: {out['matrix']} shape={pa.shape}")

        print("📊 Building per-element clade counts + stats...")

        clade_sizes = meta_df["Lineage"].value_counts().to_dict()
        if scheme == "Split_Clade1A1B" and mode == "With_Unassigned":
            expected_sizes = {"Clade 1A":169, "Clade 1B":759, "Clade 2":894,
                              "Clade 3":1357, "Unassigned":127}
            if clade_sizes != expected_sizes:
                raise SystemExit(f"[FAIL] Clade membership differs from Table S5: {clade_sizes} != {expected_sizes}")
            print("[PASS] Table S5 denominator: 3306 = 169 + 759 + 894 + 1357 + 127")

        def clade_header(c):
            return f"{c} (n={clade_sizes.get(c, 0)})"

        counts_by_clade = {}
        for c in clade_order:
            samples_c = meta_df.loc[meta_df["Lineage"] == c, "Sample"].values
            counts_by_clade[c] = pa.loc[samples_c].sum(axis=0) if len(samples_c) else pd.Series(0, index=pa.columns)

        total_present = pa.sum(axis=0)

        omnibus_rows = []
        used_clades_by_element = {}

        for symbol in pa.columns:
            pres = {c: int(counts_by_clade[c].loc[symbol]) for c in clade_order}
            ns   = {c: int(clade_sizes.get(c, 0)) for c in clade_order}
            absn = {c: ns[c] - pres[c] for c in clade_order}
            tot_pres = int(total_present.loc[symbol])

            # A clade that is 0% or 100% for a gene is CRITICAL evidence
            # of clade enrichment, not a group to remove from the hypothesis.
            tested_clades = [c for c in clade_order if ns[c] > 0]
            used_clades_by_element[symbol] = tested_clades

            pval = np.nan
            method = ""
            chi2 = np.nan
            dof = np.nan
            reason = ""

            if tot_pres < MIN_TOTAL_PRESENT:
                reason = f"rare (<{MIN_TOTAL_PRESENT} total present)"
            elif tot_pres == sum(ns.values()):
                reason = "invariant in all analyzed genomes"
            elif len(tested_clades) < 2:
                reason = "fewer than two populated clades"
            else:
                table = np.array([[pres[c], absn[c]] for c in tested_clades], dtype=int)
                try:
                    pval, method, chi2, dof = chi2_with_optional_montecarlo(table)
                except Exception as e:
                    reason = f"chi2 failed: {e}"

            omnibus_rows.append({
                "Element symbol": symbol,
                "Total_Present": tot_pres,
                "Tested_Clades": ",".join(tested_clades),
                "Test_Method": method,
                "Chi2": chi2,
                "DoF": dof,
                "Stats p-value": pval,
                "Stats FDR": np.nan,
                "Reason_Not_Tested": reason,
            })

        omnibus_df = pd.DataFrame(omnibus_rows)

        p_num = pd.to_numeric(omnibus_df["Stats p-value"], errors="coerce")
        tested_mask = p_num.notna()
        if tested_mask.sum() > 0:
            _, fdr_vals, _, _ = multipletests(p_num.loc[tested_mask].values.astype(float), method="fdr_bh")
            omnibus_df.loc[tested_mask, "Stats FDR"] = fdr_vals

        omnibus_df.to_csv(out["omnibus"], sep="\t", index=False)
        print(f"[OK] Omnibus stats: {out['omnibus']}")

        print("📊 Running post-hoc pairwise Fisher (Holm within element) for significant elements...")

        fdr_num = pd.to_numeric(omnibus_df["Stats FDR"], errors="coerce")
        sig_elements = omnibus_df.loc[tested_mask & (fdr_num <= ALPHA), "Element symbol"].tolist()

        pairwise_rows = []
        cld_letters_by_element = {e: {c: "" for c in clade_order} for e in pa.columns}

        for symbol in sig_elements:
            pres = {c: int(counts_by_clade[c].loc[symbol]) for c in clade_order}
            ns   = {c: int(clade_sizes.get(c, 0)) for c in clade_order}
            absn = {c: ns[c] - pres[c] for c in clade_order}

            tested_clades = used_clades_by_element.get(symbol, [])
            if len(tested_clades) < 2:
                continue

            pairs = list(itertools.combinations(tested_clades, 2))
            pvals = []
            tmp_rows = []

            for a, b in pairs:
                a1, a0 = pres[a], absn[a]
                b1, b0 = pres[b], absn[b]

                _, p_raw = fisher_exact([[a1, a0], [b1, b0]], alternative="two-sided")
                or05 = odds_ratio_with_haldene(a1, a0, b1, b0, add=0.5)

                tmp_rows.append({
                    "Element symbol": symbol,
                    "Group1": a,
                    "Group2": b,
                    "G1_present": a1,
                    "G1_absent": a0,
                    "G2_present": b1,
                    "G2_absent": b0,
                    "OR_0.5": or05,
                    "p_raw": float(p_raw),
                    "p_holm": np.nan,
                })
                pvals.append(float(p_raw))

            if pvals:
                _, p_holm, _, _ = multipletests(pvals, method="holm")
                for i in range(len(tmp_rows)):
                    tmp_rows[i]["p_holm"] = float(p_holm[i])

            pairwise_rows.extend(tmp_rows)

            p_adj_map = {}
            for r in tmp_rows:
                key = tuple(sorted((r["Group1"], r["Group2"])))
                p_adj_map[key] = r["p_holm"]

            cld_all = make_cld_letters(tested_clades, p_adj_map, alpha=ALPHA)
            for c in tested_clades:
                cld_letters_by_element[symbol][c] = cld_all.get(c, "")

        pairwise_df = pd.DataFrame(pairwise_rows)
        pairwise_df.to_csv(out["pairwise"], sep="\t", index=False)
        print(f"[OK] Pairwise stats: {out['pairwise']}")

        summary_rows = []
        p_map = dict(zip(omnibus_df["Element symbol"], omnibus_df["Stats p-value"]))
        f_map = dict(zip(omnibus_df["Element symbol"], omnibus_df["Stats FDR"]))

        for symbol in pa.columns:
            t = sym_to_type.get(symbol, "NA")
            nm = sym_to_name.get(symbol, "NA")

            letters = cld_letters_by_element.get(symbol, {c: "" for c in clade_order})

            row = {
                "Type": t,
                "Element symbol": symbol,
                "Element name": nm,
                "Total (n)": int(total_present.loc[symbol]),
            }
            for c in clade_order:
                n = int(clade_sizes.get(c, 0))
                cnt = int(counts_by_clade[c].loc[symbol]) if n > 0 else 0
                row[clade_header(c)] = format_cell(cnt, n, letters.get(c, ""))

            pval = p_map.get(symbol, np.nan)
            fdr  = f_map.get(symbol, np.nan)
            row["Stats p-value"] = "" if pd.isna(pval) else float(pval)
            row["Stats FDR"]     = "" if pd.isna(fdr) else float(fdr)

            summary_rows.append(row)

        summary_df = pd.DataFrame(summary_rows)

        _tested = pd.to_numeric(summary_df["Stats p-value"], errors="coerce").notna()
        _total  = pd.to_numeric(summary_df["Total (n)"], errors="coerce").fillna(0).astype(int)

        summary_df["_tested_sort"] = _tested
        summary_df["_total_sort"] = _total

        summary_df = (
            summary_df.sort_values(["_tested_sort", "_total_sort", "Type", "Element symbol"],
                                   ascending=[False, False, True, True])
                     .drop(columns=["_tested_sort", "_total_sort"])
                     .reset_index(drop=True)
        )

        summary_df.to_csv(out["summary"], sep="\t", index=False)
        summary_df.to_csv(out["table"], sep="\t", index=False)
        print(f"[OK] Summary (paper-style): {out['summary']}")
        print(f"[OK] ElementCounts TSV (same content): {out['table']}")

        # Audit original manuscript claims, without silently forcing expected
        # counts to match the current input data or overwriting old output.
        if scheme == "Split_Clade1A1B" and mode == "With_Unassigned":
            expected_s5 = {
                "sseK3": [10,70,0,0,13],
                "gogB": [0,0,0,0,22],
                "ybtP": [0,0,0,0,10],
                "gyrA_S83F": [0,3,6,0,2],
                "mdsB": [169,759,894,1356,126],
                "golT": [169,759,894,1357,118],
                "sodC1": [169,759,894,1357,121],
                "vasL": [169,759,894,1357,121],
                "gtgA": [169,759,894,1357,102],
                "pagO": [169,759,890,1357,125],
                "allD": [169,759,891,1357,127],
            }
            table_s5_audit=[]
            omnibus_lookup=omnibus_df.set_index("Element symbol")
            for gene,expected in expected_s5.items():
                observed=[int(counts_by_clade[c].get(gene,0)) for c in clade_order]
                row=omnibus_lookup.loc[gene] if gene in omnibus_lookup.index else None
                table_s5_audit.append({
                    "Gene":gene, "Expected_clade_counts":",".join(map(str,expected)),
                    "Observed_clade_counts":",".join(map(str,observed)),
                    "Counts_match":observed==expected,
                    "Omnibus_p":float(row["Stats p-value"]) if row is not None and pd.notna(row["Stats p-value"]) else np.nan,
                    "FDR_BH":float(row["Stats FDR"]) if row is not None and pd.notna(row["Stats FDR"]) else np.nan,
                    "Not_tested_reason":row["Reason_Not_Tested"] if row is not None else "not in AMRFinder data",
                    "Compact_letters":",".join(cld_letters_by_element.get(gene,{}).get(c,"") for c in clade_order),
                })
            audit=pd.DataFrame(table_s5_audit)
            audit_file=os.path.join(output_dir,"Table_S5_Flagged_Genes_Audit.tsv")
            audit.to_csv(audit_file,sep="\t",index=False)
            print("[AUDIT] Table S5 comparison:",audit_file)
            print(audit.to_string(index=False))
            if not audit["Counts_match"].all():
                print("[REVIEW REQUIRED] Some master-derived gene counts differ from pasted Table S5; see audit before manuscript edits.")

        if args.stats_only:
            # Statistics and all sample-level binary calls are retained; only
            # the huge 600-dpi figures are skipped for fast Table S5 auditing.
            matrix_with_meta = pa.copy()
            matrix_with_meta.insert(0, "Lineage", meta_df.set_index("Sample")
                                    .loc[matrix_with_meta.index, "Lineage"].astype(str).values)
            with pd.ExcelWriter(out["excel"],engine="openpyxl") as xw:
                summary_df.to_excel(xw,sheet_name="Summary_ForPaper",index=False)
                omnibus_df.to_excel(xw,sheet_name="Omnibus_Stats",index=False)
                pairwise_df.to_excel(xw,sheet_name="Pairwise_Fisher_Holm",index=False)
                meta_df.to_excel(xw,sheet_name="Sample_Clade_Metadata",index=False)
                matrix_with_meta.to_excel(xw,sheet_name="PresenceAbsence_Matrix",index=True)
                if scheme == "Split_Clade1A1B" and mode == "With_Unassigned":
                    audit.to_excel(xw,sheet_name="Table_S5_Flagged_Audit",index=False)
            print(f"[PASS] Stats-only Table S5 results: {output_dir}")
            print(f"[OUTPUT] {out['excel']}")
            continue

        print("🎨 Building heatmap (Types grouped together)...")

        type_by_sym = {sym: str(sym_to_type.get(sym, "NA")) for sym in pa.columns}
        type_groups = {}
        for sym, t in type_by_sym.items():
            type_groups.setdefault(t, []).append(sym)

        ordered_types = sorted(type_groups.keys())

        ordered_cols = []
        for t in ordered_types:
            syms = type_groups[t]
            block = pa[syms]
            ordered_cols.extend(cluster_columns_within_group(block))

        plot_data = pa[ordered_cols].copy()

        row_colors = meta_df.set_index("Sample")["Lineage"].map(clade_palette)
        type_palette = safe_type_palette(ordered_types)
        col_colors = pd.Series([type_palette[type_by_sym[sym]] for sym in ordered_cols], index=ordered_cols)

        g = sns.clustermap(
            plot_data,
            cmap=my_cmap,
            row_cluster=False,
            col_cluster=False,
            row_colors=row_colors,
            col_colors=col_colors,
            figsize=(34, 80),
            yticklabels=True,
            xticklabels=True,
            linewidths=0,
            cbar_pos=None
        )

        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=6)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=2)

        presence_legend = [
            Patch(facecolor="#FF0000", label="Present (1)"),
            Patch(facecolor="#FFFFFF", edgecolor="gray", label="Absent (0)"),
        ]
        clade_legend = [Patch(facecolor=clade_palette[c], label=c) for c in clade_order if clade_sizes.get(c, 0) > 0]
        type_legend  = [Patch(facecolor=type_palette[t], label=t) for t in ordered_types]

        leg1 = plt.legend(handles=presence_legend, title="Presence", loc="upper left", bbox_to_anchor=(1.05, 1.00))
        plt.gca().add_artist(leg1)
        leg2 = plt.legend(handles=clade_legend, title="Clade", loc="upper left", bbox_to_anchor=(1.05, 0.86))
        plt.gca().add_artist(leg2)
        plt.legend(handles=type_legend, title="Type", loc="upper left", bbox_to_anchor=(1.05, 0.60))

        g.savefig(out["heat_pdf"], bbox_inches="tight")
        g.savefig(out["heat_png"], dpi=600, bbox_inches="tight")
        plt.close("all")

        print(f"[OK] Heatmap PDF: {out['heat_pdf']}")
        print(f"[OK] Heatmap PNG: {out['heat_png']}")

        print("💾 Writing Excel workbook...")

        meta_out = meta_df.copy()
        meta_out["Lineage"] = meta_out["Lineage"].astype(str)

        matrix_with_meta = pa.copy()
        matrix_with_meta.insert(
            0,
            "Lineage",
            meta_out.set_index("Sample").loc[matrix_with_meta.index, "Lineage"].values
        )

        with pd.ExcelWriter(out["excel"], engine="openpyxl") as xw:
            summary_df.to_excel(xw, sheet_name="Summary_ForPaper", index=False)
            omnibus_df.to_excel(xw, sheet_name="Omnibus_Stats", index=False)
            pairwise_df.to_excel(xw, sheet_name="Pairwise_Fisher_Holm", index=False)
            meta_out.to_excel(xw, sheet_name="Sample_Clade_Metadata", index=False)
            matrix_with_meta.to_excel(xw, sheet_name="PresenceAbsence_Matrix", index=True)

        if os.path.exists(tmp_badlines):
            try:
                import shutil
                shutil.copy2(tmp_badlines, out["badlines"])
            except Exception:
                pass

        print("✅ DONE")
        print(f"[OK] Excel: {out['excel']}")
        print("\n[INFO] Quick counts:")
        print(meta_df["Lineage"].value_counts().reindex(clade_order).fillna(0).astype(int))
        print(f"[INFO] Significant elements (FDR <= {ALPHA}): {len(sig_elements)}")

print("\n🎉 ALL DONE")
print(f"[OK] Results root: {ROOT_OUT_DIR}")
