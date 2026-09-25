#!/usr/bin/env python3
import subprocess
import sys
import os
import re
from pathlib import Path
import itertools
import warnings

# =============================================================================
# 1) DEPENDENCY MANAGER
# =============================================================================
def check_and_install_dependencies():
    required = {
        "pandas": "pandas",
        "numpy": "numpy",
        "matplotlib": "matplotlib",
        "scipy": "scipy",
        "statsmodels": "statsmodels",
        "openpyxl": "openpyxl",
    }
    for imp, pip_name in required.items():
        try:
            __import__(imp)
        except ImportError:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])

check_and_install_dependencies()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
# CONFIG
# =============================================================================
REPO_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = str(Path(os.environ.get("SEN_ROOT", REPO_ROOT)))
SEN_GENOMES = os.environ.get("SEN_METADATA_FILE", os.path.join(BASE_DIR, "metadata", "SEN_Genomes.csv"))
COLLECTION_YEAR_ITOL = os.path.join(BASE_DIR, "itol_4_collection_year.txt")

metadata_path_candidates = [
    os.environ.get("SEN_CLADE_METADATA", os.path.join(BASE_DIR, "metadata", "final_clade_metadata.tsv")),
    os.path.join(BASE_DIR, "iqtree_final", "ITOL", "Clade_metadata.txt"),
    os.path.join(BASE_DIR, "ITOL", "Clade_metadata.txt"),
    os.path.join(BASE_DIR, "iqtree", "ITOL", "Clade_metadata.txt"),
]
CLADE_META = next((p for p in metadata_path_candidates if os.path.exists(p)), None)
if CLADE_META is None:
    raise SystemExit("[ERROR] Could not find Clade_metadata.txt in expected locations.")

if not os.path.exists(COLLECTION_YEAR_ITOL):
    raise SystemExit(f"[ERROR] Could not find collection year annotation file: {COLLECTION_YEAR_ITOL}")

ROOT_OUT_DIR = os.path.join(BASE_DIR, "iqtree_final", "clade_summary")
os.makedirs(ROOT_OUT_DIR, exist_ok=True)

# Stats
ALPHA = 0.05
EXCLUDE_ZERO_OR_ALL = True
MIN_TOTAL_CATEGORY = 5

# Plot styling
YMAX = 100
TITLE_FONTSIZE = 20
AXISLABEL_FONTSIZE = 14
TICK_FONTSIZE = 12
SIGNIF_BASE_FONTSIZE = 12
SIGNIF_FONTSIZE = SIGNIF_BASE_FONTSIZE * 1.25
SIGNIF_FONTWEIGHT = "bold"

# -----------------------------------------------------------------------------
# CATEGORY COLOR MAPS
# -----------------------------------------------------------------------------
SOURCE_COLORS = {
    "Chicken-US": "#d62728",
    "Human": "#2ca02c",
    "International": "#17becf",
    "Other-US": "#9467bd",
    "Other": "#7f7f7f",
    "Unknown": "#c7c7c7",
}

QR_COLORS = {
    "QS": "#d62728",
    "S83Y": "#1f77b4",
    "D87Y": "#9467bd",
    "D87N": "#2ca02c",
    "D87G": "#ff7f0e",
    "Others": "#7f7f7f",
    "Other": "#7f7f7f",
    "Unknown": "#c7c7c7",
}

YEAR_COLORS = {
    "<=2000": "#e6194b",
    "2001-2005": "#3cb44b",
    "2006-2010": "#ffe119",
    "2011-2015": "#4363d8",
    "2016-2020": "#f58231",
    "2021-2025": "#911eb4",
    "Unknown": "#c7c7c7",
}

LEGEND_TITLES = {
    "Clean_Source": "Sources",
    "QR_marker": "QR markers",
    "Collection_Year": "Collection year",
}

# =============================================================================
# HELPERS
# =============================================================================
def normalize_sample_id(x: str) -> str:
    if pd.isna(x):
        return x
    s = str(x).strip()
    s = os.path.basename(s)
    s = re.sub(r"\.tsv$|\.txt$|\.csv$", "", s, flags=re.IGNORECASE)
    s = s.strip("_")
    s = re.sub(r"(_trimmed.*)$", "", s)
    s = re.sub(r"(_pure.*)$", "", s)
    m = re.search(r"(SRR\d+)", s)
    return m.group(1) if m else s

def robust_read_csv(path, sep=","):
    try:
        return pd.read_csv(path, sep=sep, engine="c", low_memory=False)
    except Exception:
        return pd.read_csv(path, sep=sep, engine="python", on_bad_lines="skip")

def save_table_both(df: pd.DataFrame, base_path_no_ext: str, index: bool = False):
    df.to_csv(base_path_no_ext + ".tsv", sep="\t", index=index)
    df.to_csv(base_path_no_ext + ".csv", index=index)

def clean_source(x: str) -> str:
    if pd.isna(x) or str(x).strip() == "":
        return "Unknown"
    s = str(x).strip()
    if "NON-US" in s.upper():
        return "International"
    return s

def clean_qr_marker(qr_group, amr_pattern) -> str:
    q = "" if pd.isna(qr_group) else str(qr_group).strip().upper()
    if ("NON" in q) or (q in {"0", "NO", "FALSE", "N"}):
        return "QS"
    raw = "" if pd.isna(amr_pattern) else str(amr_pattern).strip().upper()
    if raw in {"", "NAN", "NA", "NONE"}:
        return "QS"
    clean = (
        raw.replace("GYRA_", "").replace("GYRA", "")
           .replace("PARC_", "").replace("PARC", "")
           .replace("PARG_", "").replace("PARG", "")
           .replace(" ", "")
    )
    primary = {"S83Y", "D87Y", "D87N", "D87G"}
    tokens = re.split(r"[;|,+/]+", clean)
    tokens = [t for t in tokens if t]
    for t in tokens:
        if t in primary:
            return t
    return "Others" if tokens else "QS"

def parse_collection_year_itol(path: str) -> pd.DataFrame:
    rows = []
    data_started = False

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.upper() == "DATA":
                data_started = True
                continue
            if not data_started:
                continue

            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 3:
                continue

            sample = normalize_sample_id(parts[0])
            color = parts[1]
            label = parts[2]

            if sample:
                rows.append({
                    "Sample": sample,
                    "Collection_Year": label,
                    "Collection_Year_Color": color
                })

    out = pd.DataFrame(rows).drop_duplicates(subset=["Sample"], keep="first")
    return out

def group_rare_categories(series: pd.Series, min_total: int = 5, other_label: str = "Other") -> pd.Series:
    vc = series.value_counts(dropna=False)
    keep = set(vc[vc >= min_total].index.tolist())
    return series.apply(lambda x: x if x in keep else other_label)

def chi2_with_optional_montecarlo(table: np.ndarray):
    chi2, p, dof, expected = chi2_contingency(table, correction=False)
    method = "chi-square"
    if np.any(expected < 5):
        try:
            chi2_mc, p_mc, dof_mc, _ = chi2_contingency(
                table, correction=False, method="montecarlo", num_resamples=5000
            )
            return float(p_mc), "chi-square (montecarlo)", float(chi2_mc), int(dof_mc)
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
        return pd.notna(p) and p <= alpha

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
            remaining.discard(g)

    return {g: "".join(sorted(assign[g])) for g in groups}

def category_order_for_attr(attr_col: str, present_cats):
    present = set(present_cats)

    if attr_col == "QR_marker":
        preferred = ["QS", "S83Y", "D87Y", "D87N", "D87G", "Others", "Other", "Unknown"]
        return [c for c in preferred if c in present] + sorted([c for c in present if c not in set(preferred)])

    if attr_col == "Clean_Source":
        preferred = ["Chicken-US", "Human", "International", "Other-US", "Other", "Unknown"]
        return [c for c in preferred if c in present] + sorted([c for c in present if c not in set(preferred)])

    if attr_col == "Collection_Year":
        preferred = ["<=2000", "2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025", "Unknown"]
        return [c for c in preferred if c in present] + sorted([c for c in present if c not in set(preferred)])

    return sorted(list(present))

def get_color_map(attr_col: str):
    if attr_col == "QR_marker":
        return QR_COLORS
    if attr_col == "Collection_Year":
        return YEAR_COLORS
    return SOURCE_COLORS

def get_legend_title(attr_col: str):
    return LEGEND_TITLES.get(attr_col, attr_col)

def safe_div(n, d):
    return (100.0 * n / d) if d and d > 0 else np.nan

def recode_lineages(df: pd.DataFrame, scheme: str) -> pd.DataFrame:
    out = df.copy()

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

def get_focus_clades(scheme: str, include_unassigned: bool):
    if scheme == "Combined_Clade1":
        base = ["Clade 1", "Clade 2", "Clade 3"]
    elif scheme == "Split_Clade1A1B":
        base = ["Clade 1A", "Clade 1B", "Clade 2", "Clade 3"]
    else:
        raise ValueError(f"Unknown scheme: {scheme}")

    return base + (["Unassigned"] if include_unassigned else [])

# =============================================================================
# ENRICHMENT ANALYSIS + PLOT
# =============================================================================
def run_attribute_enrichment(df_all, clade, attr_col, title_prefix, out_prefix, sub_dirs):
    df = df_all.copy()
    df = df[df[attr_col].notna()].copy()

    totals = df[attr_col].value_counts().to_dict()
    df_in = df[df["Lineage"] == clade].copy()
    in_counts = df_in[attr_col].value_counts().to_dict()

    categories = category_order_for_attr(attr_col, totals.keys())
    n_in_clade = int(len(df_in))

    rows = []
    for cat in categories:
        denom = int(totals.get(cat, 0))
        num = int(in_counts.get(cat, 0))
        share = safe_div(num, denom)
        within = safe_div(num, n_in_clade)
        rows.append({
            "Clade": clade,
            "Attribute": attr_col,
            "Category": cat,
            "In_Clade_n": num,
            "Total_Category_n": denom,
            "Within_Clade_%": within,
            "Share_of_Category_in_Clade_%": share,
            "Letters": "",
            "Tested_in_stats": ""
        })

    out_table = pd.DataFrame(rows)

    eligible = [
        cat for cat in categories
        if int(totals.get(cat, 0)) > 0 and (
            (not EXCLUDE_ZERO_OR_ALL) or (0 < int(in_counts.get(cat, 0)) < int(totals.get(cat, 0)))
        )
    ]

    omnibus = {
        "Clade": clade,
        "Attribute": attr_col,
        "N_in_clade": n_in_clade,
        "N_total_nonNA": int(len(df)),
        "N_categories_total": len(categories),
        "N_categories_tested": len(eligible),
        "Method": "",
        "Chi2": np.nan,
        "DoF": np.nan,
        "p_value": np.nan,
        "Reason_Not_Tested": ""
    }

    if eligible:
        tmp_sub = out_table[out_table["Category"].isin(eligible)].copy()
        omnibus["Within_Clade_%_min"] = float(np.nanmin(tmp_sub["Within_Clade_%"].values))
        omnibus["Within_Clade_%_max"] = float(np.nanmax(tmp_sub["Within_Clade_%"].values))
        omnibus["Share_of_Category_in_Clade_%_min"] = float(np.nanmin(tmp_sub["Share_of_Category_in_Clade_%"].values))
        omnibus["Share_of_Category_in_Clade_%_max"] = float(np.nanmax(tmp_sub["Share_of_Category_in_Clade_%"].values))
    else:
        omnibus["Within_Clade_%_min"] = np.nan
        omnibus["Within_Clade_%_max"] = np.nan
        omnibus["Share_of_Category_in_Clade_%_min"] = np.nan
        omnibus["Share_of_Category_in_Clade_%_max"] = np.nan

    pairwise_rows = []
    letters_map = {cat: "" for cat in categories}

    if len(eligible) >= 2:
        table = np.array([
            [int(in_counts.get(cat, 0)), int(totals.get(cat, 0)) - int(in_counts.get(cat, 0))]
            for cat in eligible
        ], dtype=int)

        try:
            pval, method, chi2, dof = chi2_with_optional_montecarlo(table)
            omnibus.update({"Method": method, "Chi2": chi2, "DoF": dof, "p_value": pval})

            pairs = list(itertools.combinations(eligible, 2))
            pvals = []
            tmp = []

            for a, b in pairs:
                a1 = int(in_counts.get(a, 0))
                a0 = int(totals.get(a, 0)) - a1
                b1 = int(in_counts.get(b, 0))
                b0 = int(totals.get(b, 0)) - b1

                _, p_raw = fisher_exact([[a1, a0], [b1, b0]])

                a_share = safe_div(a1, (a1 + a0))
                b_share = safe_div(b1, (b1 + b0))
                a_within = safe_div(a1, n_in_clade)
                b_within = safe_div(b1, n_in_clade)

                or_h = odds_ratio_with_haldene(a1, a0, b1, b0, add=0.5)
                direction = "Cat1>Cat2" if (pd.notna(or_h) and or_h > 1) else ("Cat2>Cat1" if (pd.notna(or_h) and or_h < 1) else "Equal")

                tmp.append({
                    "Clade": clade,
                    "Attribute": attr_col,
                    "Category1": a,
                    "Category2": b,
                    "Cat1_Total_n": int(totals.get(a, 0)),
                    "Cat2_Total_n": int(totals.get(b, 0)),
                    "Cat1_In_Clade_n": a1,
                    "Cat2_In_Clade_n": b1,
                    "Cat1_Share_of_Category_in_Clade_%": a_share,
                    "Cat2_Share_of_Category_in_Clade_%": b_share,
                    "Diff_Share_pp": (a_share - b_share) if (pd.notna(a_share) and pd.notna(b_share)) else np.nan,
                    "Cat1_Within_Clade_%": a_within,
                    "Cat2_Within_Clade_%": b_within,
                    "Diff_Within_pp": (a_within - b_within) if (pd.notna(a_within) and pd.notna(b_within)) else np.nan,
                    "OR_Haldane": float(or_h) if pd.notna(or_h) else np.nan,
                    "Direction": direction,
                    "p_raw": float(p_raw),
                })
                pvals.append(float(p_raw))

            if pvals:
                _, p_holm, _, _ = multipletests(pvals, method="holm")

                p_adj_map = {}
                for i in range(len(tmp)):
                    tmp[i]["p_holm"] = float(p_holm[i])
                    p_adj_map[tuple(sorted((tmp[i]["Category1"], tmp[i]["Category2"])))] = float(p_holm[i])

                pairwise_rows = tmp

                cld = make_cld_letters(eligible, p_adj_map, alpha=ALPHA)
                for cat in eligible:
                    letters_map[cat] = cld.get(cat, "")

        except Exception as e:
            omnibus["Reason_Not_Tested"] = str(e)
    else:
        omnibus["Reason_Not_Tested"] = "Fewer than 2 eligible categories for testing."

    out_table["Letters"] = out_table["Category"].map(letters_map).fillna("")
    out_table["Tested_in_stats"] = out_table["Category"].apply(lambda x: "YES" if x in set(eligible) else "")
    out_table = out_table.sort_values("Share_of_Category_in_Clade_%", ascending=False).reset_index(drop=True)

    tag = clade.replace(" ", "")
    save_table_both(out_table, os.path.join(sub_dirs["ENR_TABLES"], f"{out_prefix}_{tag}_table"))
    save_table_both(pd.DataFrame([omnibus]), os.path.join(sub_dirs["ENR_TABLES"], f"{out_prefix}_{tag}_omnibus"))
    save_table_both(pd.DataFrame(pairwise_rows), os.path.join(sub_dirs["ENR_TABLES"], f"{out_prefix}_{tag}_pairwise"))

    color_map = get_color_map(attr_col)

    plt.figure(figsize=(max(10, 0.75 * len(out_table)), 6))
    x_pos = np.arange(len(out_table))
    y_vals = out_table["Share_of_Category_in_Clade_%"].astype(float).values

    plt.bar(
        x_pos, y_vals,
        color=[color_map.get(c, "#bdbdbd") for c in out_table["Category"]],
        edgecolor="white", linewidth=0.5
    )

    plt.xticks(
        x_pos,
        [f"{c} (n={int(n)})" for c, n in zip(out_table["Category"], out_table["Total_Category_n"])],
        rotation=45, ha="right", fontsize=TICK_FONTSIZE
    )
    plt.ylabel("Share of category in clade (%)", fontsize=AXISLABEL_FONTSIZE)
    plt.title(f"{clade} – {title_prefix} enrichment", fontsize=TITLE_FONTSIZE)
    plt.ylim(0, YMAX)

    pad = max(0.01 * (np.nanmax(y_vals) if len(y_vals) > 0 else 1), 0.6)
    for i, (val, let) in enumerate(zip(y_vals, out_table["Letters"])):
        if let:
            plt.text(i, min(val + pad, YMAX - 1), let, ha="center", va="bottom",
                     fontsize=SIGNIF_FONTSIZE, fontweight=SIGNIF_FONTWEIGHT)

    plt.legend(
        handles=[mpl.patches.Patch(color=color_map.get(c, "#bdbdbd"), label=c) for c in categories],
        title=get_legend_title(attr_col),
        bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False
    )

    plt.tight_layout()
    out_png = os.path.join(sub_dirs["ENR_CHARTS"], f"{out_prefix}_{tag}.png")
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.savefig(out_png.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()

    return out_table, pd.DataFrame([omnibus]), pd.DataFrame(pairwise_rows), out_png

# =============================================================================
# COMPOSITION STACKED BARS
# =============================================================================
def make_composition_stacked(df_all, attr_col, title, out_prefix, focus_clades, sub_dirs):
    df = df_all[df_all["Lineage"].isin(focus_clades) & df_all[attr_col].notna()].copy()
    if df.empty:
        return None, None, None

    color_map = get_color_map(attr_col)
    legend_title = get_legend_title(attr_col)
    cats = category_order_for_attr(attr_col, df[attr_col].unique())

    clade_sizes = df["Lineage"].value_counts().to_dict()
    total_cat_global = df[attr_col].value_counts().to_dict()

    rows = []
    for cl in focus_clades:
        denom_clade = int(clade_sizes.get(cl, 0))
        vc = df[df["Lineage"] == cl][attr_col].value_counts().to_dict()

        for cat in cats:
            n = int(vc.get(cat, 0))
            within_pct = safe_div(n, denom_clade)
            denom_cat = int(total_cat_global.get(cat, 0))
            share_cat_in_clade = safe_div(n, denom_cat)

            rows.append({
                "Clade": cl,
                "Attribute": attr_col,
                "Category": cat,
                "In_Clade_n": n,
                "Within_Clade_%": within_pct,
                "Share_of_Category_in_Clade_%": share_cat_in_clade,
                "Total_Category_n": denom_cat,
                "Clade_n": denom_clade
            })

    comp_df = pd.DataFrame(rows)
    comp_base = os.path.join(sub_dirs["COMP_TABLES"], f"{out_prefix}_composition_table")
    save_table_both(comp_df, comp_base)

    plt.figure(figsize=(10, 6))
    x_pos = np.arange(len(focus_clades))
    bottom = np.zeros(len(focus_clades))

    for cat in cats:
        vals = []
        for cl in focus_clades:
            v = comp_df[(comp_df["Clade"] == cl) & (comp_df["Category"] == cat)]["Within_Clade_%"].values
            vals.append(float(v[0]) if len(v) else 0.0)

        plt.bar(
            x_pos, vals,
            bottom=bottom,
            color=color_map.get(cat, "#bdbdbd"),
            label=cat,
            edgecolor="white", linewidth=0.4
        )
        bottom += np.array(vals)

    plt.xticks(
        x_pos,
        [f"{cl}\n(n={int(clade_sizes.get(cl, 0))})" for cl in focus_clades],
        fontsize=TICK_FONTSIZE
    )
    plt.ylabel("Within-clade composition (%)", fontsize=AXISLABEL_FONTSIZE)
    plt.title(title, fontsize=TITLE_FONTSIZE)
    plt.ylim(0, 100)

    plt.legend(
        handles=[mpl.patches.Patch(color=color_map.get(c, "#bdbdbd"), label=c) for c in cats],
        title=legend_title,
        bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False
    )

    plt.tight_layout()
    out_png = os.path.join(sub_dirs["COMP_CHARTS"], f"{out_prefix}_composition.png")
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.savefig(out_png.replace(".png", ".pdf"), bbox_inches="tight")
    plt.close()

    return comp_df, comp_base + ".tsv", out_png

# =============================================================================
# MAIN PIPELINE
# =============================================================================
print("[1] Initializing Data...")
gen = robust_read_csv(SEN_GENOMES, sep=",")
gen.columns = [c.strip() for c in gen.columns]

clade_df = robust_read_csv(CLADE_META, sep="\t")
clade_df["Sample"] = clade_df["Sample"].map(normalize_sample_id)
clade_df = clade_df[clade_df["Lineage"] != "Reference"].copy()

year_df = parse_collection_year_itol(COLLECTION_YEAR_ITOL)

run_col = next((c for c in ["Run", "Sample", "Accession", "SRR"] if c in gen.columns), None)
if not run_col:
    raise SystemExit("[ERROR] SEN_Genomes.csv missing Sample ID column.")
gen["Sample"] = gen[run_col].map(normalize_sample_id)

merged_raw = clade_df.merge(gen, on="Sample", how="left")
merged_raw = merged_raw.merge(year_df[["Sample", "Collection_Year"]], on="Sample", how="left")

merged_raw["Clean_Source"] = merged_raw["Source"].map(clean_source) if "Source" in merged_raw.columns else "Unknown"
merged_raw["QR_marker"] = merged_raw.apply(
    lambda r: clean_qr_marker(r.get("QR_group", np.nan), r.get("AMR_pattern", np.nan)),
    axis=1
)
merged_raw["Collection_Year"] = merged_raw["Collection_Year"].fillna("Unknown")

merged_raw["Clean_Source"] = group_rare_categories(merged_raw["Clean_Source"], MIN_TOTAL_CATEGORY, "Other")
merged_raw["QR_marker"] = group_rare_categories(merged_raw["QR_marker"], MIN_TOTAL_CATEGORY, "Others")
valid_year_bins = {"<=2000", "2001-2005", "2006-2010", "2011-2015", "2016-2020", "2021-2025", "Unknown"}
merged_raw["Collection_Year"] = merged_raw["Collection_Year"].apply(lambda x: x if x in valid_year_bins else "Unknown")

print("\n[DEBUG] Unique lineage labels found in metadata:")
print(sorted(merged_raw["Lineage"].dropna().unique()))

print("\n[DEBUG] Unique collection year bins found:")
print(sorted(merged_raw["Collection_Year"].dropna().unique()))

schemes = ["Combined_Clade1", "Split_Clade1A1B"]
unassigned_modes = ["With_Unassigned", "Without_Unassigned"]

for scheme in schemes:
    print(f"\n==============================")
    print(f"🚀 Processing scheme: {scheme}")
    print(f"==============================")

    scheme_df = recode_lineages(merged_raw, scheme)

    print("[DEBUG] Unique lineage labels after recoding:")
    print(sorted(scheme_df["Lineage"].dropna().unique()))

    for mode in unassigned_modes:
        print(f"\n➡ Processing mode: {mode}")
        folder_path = os.path.join(ROOT_OUT_DIR, scheme, mode)

        sub_dirs = {
            "ENR_TABLES": os.path.join(folder_path, "enrichment", "Tables"),
            "ENR_CHARTS": os.path.join(folder_path, "enrichment", "Charts"),
            "COMP_TABLES": os.path.join(folder_path, "composition", "Tables"),
            "COMP_CHARTS": os.path.join(folder_path, "composition", "Charts"),
        }
        for d in sub_dirs.values():
            os.makedirs(d, exist_ok=True)

        include_unassigned = (mode == "With_Unassigned")
        current_focus = get_focus_clades(scheme, include_unassigned)

        if include_unassigned:
            current_df = scheme_df.copy()
            current_df["Lineage"] = current_df["Lineage"].fillna("Unassigned")
        else:
            current_df = scheme_df[scheme_df["Lineage"] != "Unassigned"].copy()

        current_df["Lineage"] = pd.Categorical(current_df["Lineage"], categories=current_focus, ordered=True)

        save_table_both(current_df, os.path.join(folder_path, "Merged_Metadata"))

        enrichment_results = []
        composition_results = []

        attribute_specs = [
            ("Clean_Source", "Source", "Enrichment_Source", "Composition – Sources", "Composition_Source"),
            ("QR_marker", "QR marker", "Enrichment_QR", "Composition – QR markers", "Composition_QR"),
            ("Collection_Year", "Collection year", "Enrichment_CollectionYear", "Composition – Collection year", "Composition_CollectionYear"),
        ]

        for cl in current_focus:
            for attr, label, pref_enr, _, _ in attribute_specs:
                enr_table, enr_omnibus, enr_pairwise, _ = run_attribute_enrichment(
                    current_df, cl, attr, label, pref_enr, sub_dirs
                )
                enrichment_results.append({
                    "clade": cl,
                    "attr": attr,
                    "table": enr_table,
                    "omnibus": enr_omnibus,
                    "pairwise": enr_pairwise,
                })

        for attr, _, _, comp_title, pref_comp in attribute_specs:
            comp_df, _, _ = make_composition_stacked(current_df, attr, comp_title, pref_comp, current_focus, sub_dirs)
            composition_results.append({
                "attr": attr,
                "table": comp_df
            })

        excel_path = os.path.join(folder_path, f"Summary_{scheme}_{mode}.xlsx")
        with pd.ExcelWriter(excel_path, engine="openpyxl") as xw:
            current_df.to_excel(xw, sheet_name="Merged_Metadata", index=False)

            clade_sizes_df = (
                current_df["Lineage"].value_counts()
                .reindex(current_focus).fillna(0)
                .reset_index()
            )
            clade_sizes_df.columns = ["Lineage", "N"]
            clade_sizes_df.to_excel(xw, sheet_name="CladeSizes", index=False)

            # composition sheets
            for comp in composition_results:
                attr = comp["attr"]
                comp_df = comp["table"]
                if comp_df is None or comp_df.empty:
                    continue

                if attr == "Clean_Source":
                    sheet_name = "Source_Composition"
                elif attr == "QR_marker":
                    sheet_name = "QR_Composition"
                elif attr == "Collection_Year":
                    sheet_name = "CollectionYear_Composition"
                else:
                    sheet_name = f"{attr[:25]}_Composition"

                comp_df.to_excel(xw, sheet_name=sheet_name[:31], index=False)

            # enrichment sheets
            for res in enrichment_results:
                clade = res["clade"].replace(" ", "")
                attr = res["attr"]

                if attr == "Clean_Source":
                    prefix = "Source"
                elif attr == "QR_marker":
                    prefix = "QR"
                elif attr == "Collection_Year":
                    prefix = "CollectionYear"
                else:
                    prefix = attr[:12]

                sheet_base = f"{prefix}_{clade}"

                if res["table"] is not None and not res["table"].empty:
                    res["table"].to_excel(xw, sheet_name=f"{sheet_base}_Tbl"[:31], index=False)
                if res["omnibus"] is not None and not res["omnibus"].empty:
                    res["omnibus"].to_excel(xw, sheet_name=f"{sheet_base}_Omni"[:31], index=False)
                if res["pairwise"] is not None and not res["pairwise"].empty:
                    res["pairwise"].to_excel(xw, sheet_name=f"{sheet_base}_Pair"[:31], index=False)

print("\nDONE ✅ Folders created in:", ROOT_OUT_DIR)
