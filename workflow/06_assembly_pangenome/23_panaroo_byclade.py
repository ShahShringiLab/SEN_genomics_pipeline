#!/usr/bin/env python3
import subprocess
import sys
import os
import time
import re
import warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

warnings.simplefilter("ignore")

# --- Dependency Check & Auto-Install ---
def check_and_install():
    req = ["pandas", "numpy", "scipy", "statsmodels", "pyarrow", "openpyxl", "matplotlib"]
    for r in req:
        try:
            __import__(r if r != "statsmodels" else "statsmodels")
        except ImportError:
            print(f"📦 Installing missing package: {r} ...")
            subprocess.check_call([sys.executable, "-m", "pip", "install", r])

check_and_install()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import hypergeom, fisher_exact
from statsmodels.stats.multitest import multipletests

# -------------------------
# CONFIGURATION
# -------------------------
REPO_ROOT      = Path(__file__).resolve().parents[2]
WD             = Path(os.environ.get("SEN_ROOT", REPO_ROOT))
PA_MATRIX_PATH = os.environ.get("SEN_PANAROO_MATRIX", str(WD / "Panaroo_Run" / "results" / "gene_presence_absence.csv"))
PA_DATA_PATH   = os.environ.get("SEN_PANAROO_GENE_DATA", str(WD / "Panaroo_Run" / "results" / "gene_data.csv"))
METADATA_PATH  = os.environ.get("SEN_CLADE_METADATA", str(WD / "metadata" / "final_clade_metadata.tsv"))
OUTROOT        = os.environ.get("SEN_PANAROO_CLADE_OUT", str(WD / "iqtree_final" / "Panaroo_byclade_robust"))

ACCESSORY_MIN_FREQ = 0.01
ACCESSORY_MAX_FREQ = 0.99
FDR_THRESHOLD      = 0.05
UNASSIGNED_LABEL   = "Unassigned"
REFERENCE_LABEL    = "Reference"

# Clade-definer thresholds for clade-vs-rest
CLADEDEF_IN_PCT    = 90.0
CLADEDEF_OUT_PCT   = 5.0

# Pangenome category thresholds
CORE_PCT           = 95.0
SHELL_PCT          = 15.0

# Openness / accumulation
ACCUM_REPS         = 100
RNG_SEED           = 123

EXCEL_MAX_ROWS     = 1048576

os.makedirs(OUTROOT, exist_ok=True)

# -------------------------
# HELPERS
# -------------------------
def normalize_sample_id(x: str) -> str:
    if pd.isna(x):
        return x
    s = str(x).strip()
    s = os.path.basename(s)
    s = re.sub(r"\.tsv$|\.txt$|\.csv$|\.fa$|\.fasta$|\.gff$", "", s, flags=re.IGNORECASE)
    s = s.replace("_trimmed", "").replace("_pure", "")
    s = s.strip("_")
    m = re.search(r"(SRR\d+)", s)
    return m.group(1) if m else s

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
    elif scheme == "Split_Clade1A1B":
        clade_order = ["Clade 1A", "Clade 1B", "Clade 2", "Clade 3"]
    else:
        raise ValueError(f"Unknown scheme: {scheme}")

    if include_unassigned:
        clade_order = clade_order + ["Unassigned"]

    return clade_order, include_unassigned

def load_gene_info():
    if not os.path.exists(PA_DATA_PATH):
        return pd.DataFrame(columns=["Gene", "Full_Description"])

    df = pd.read_csv(PA_DATA_PATH, usecols=["clustering_id", "description"], engine="pyarrow")
    df = df.drop_duplicates(subset=["clustering_id"])
    df.columns = ["Gene", "Full_Description"]
    return df

def find_sample_columns(pa_df: pd.DataFrame, sample_set: set):
    valid_cols = []
    col_to_sample = {}
    for c in pa_df.columns:
        nc = normalize_sample_id(c)
        if nc in sample_set:
            valid_cols.append(c)
            col_to_sample[c] = nc
    return valid_cols, col_to_sample

def save_tsv_csv(df: pd.DataFrame, prefix_no_ext: str):
    df.to_csv(prefix_no_ext + ".csv", index=False)
    df.to_csv(prefix_no_ext + ".tsv", sep="\t", index=False)

def save_excel_safe(outdir: str, sheet_map: dict):
    xlsx = os.path.join(outdir, "Panaroo_Robust_Summary.xlsx")
    with pd.ExcelWriter(xlsx, engine="openpyxl") as xw:
        for sheet_name, df in sheet_map.items():
            safe_name = sheet_name[:31]
            if len(df) <= EXCEL_MAX_ROWS:
                df.to_excel(xw, sheet_name=safe_name, index=False)
            else:
                note_df = pd.DataFrame({
                    "Message": [
                        f"Sheet '{sheet_name}' omitted because it exceeds Excel row limit.",
                        f"Rows: {len(df)}"
                    ]
                })
                note_df.to_excel(xw, sheet_name=(safe_name[:25] + "_NOTE")[:31], index=False)

def make_dir(path: str):
    os.makedirs(path, exist_ok=True)

def fit_heaps_alpha(accum_df: pd.DataFrame):
    """
    Fits log10(total_genes_mean) ~ intercept + alpha*log10(genomes_sampled)
    Returns alpha as openness proxy.
    """
    df = accum_df.copy()
    df = df[(df["Genomes_Sampled"] > 1) & (df["Total_Genes_Mean"] > 0)].copy()
    if len(df) < 2:
        return np.nan, np.nan
    x = np.log10(df["Genomes_Sampled"].astype(float).values)
    y = np.log10(df["Total_Genes_Mean"].astype(float).values)
    try:
        alpha, intercept = np.polyfit(x, y, 1)
        return float(alpha), float(intercept)
    except Exception:
        return np.nan, np.nan

def plot_pangenome_categories(summary_df: pd.DataFrame, out_png: str, out_pdf: str, title: str):
    if summary_df.empty:
        return

    df = summary_df.copy()
    x = np.arange(len(df))
    core = df["Core_Genes"].values
    shell = df["Shell_Genes"].values
    cloud = df["Cloud_Genes"].values

    plt.figure(figsize=(10, 6))
    plt.bar(x, core, label="Core")
    plt.bar(x, shell, bottom=core, label="Shell")
    plt.bar(x, cloud, bottom=core + shell, label="Cloud")
    plt.xticks(x, df["Clade"], rotation=45, ha="right")
    plt.ylabel("Number of genes")
    plt.xlabel("Clade")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

def plot_accumulation_curves(accum_df: pd.DataFrame, out_png: str, out_pdf: str, title: str):
    if accum_df.empty:
        return

    plt.figure(figsize=(9, 6))
    for clade, sub in accum_df.groupby("Clade"):
        plt.plot(sub["Genomes_Sampled"], sub["Total_Genes_Mean"], label=f"{clade} total")
    plt.xlabel("Genomes sampled")
    plt.ylabel("Mean cumulative genes")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

def plot_core_decay_curves(accum_df: pd.DataFrame, out_png: str, out_pdf: str, title: str):
    if accum_df.empty:
        return

    plt.figure(figsize=(9, 6))
    for clade, sub in accum_df.groupby("Clade"):
        plt.plot(sub["Genomes_Sampled"], sub["Core_Genes_Mean"], label=f"{clade} core")
    plt.xlabel("Genomes sampled")
    plt.ylabel("Mean core genes")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

# -------------------------
# ANALYSIS BLOCKS
# -------------------------
def build_binary_matrix_for_task(meta_df: pd.DataFrame, matrix_df: pd.DataFrame):
    valid_sample_set = set(meta_df["Sample"])
    valid_sample_cols, col_to_sample = find_sample_columns(matrix_df, valid_sample_set)

    if not valid_sample_cols:
        return None, None, None

    bin_subset = matrix_df[valid_sample_cols].notna().astype(np.uint8)
    bin_subset.columns = [col_to_sample[c] for c in valid_sample_cols]

    if bin_subset.columns.duplicated().any():
        bin_subset = bin_subset.T.groupby(level=0).max().T

    ordered_samples = meta_df["Sample"].tolist()
    bin_subset = bin_subset.reindex(columns=ordered_samples, fill_value=0)

    counts_total = bin_subset.sum(axis=1)
    M = len(ordered_samples)
    freq_global = counts_total / M if M > 0 else 0
    mask_acc = (freq_global >= ACCESSORY_MIN_FREQ) & (freq_global <= ACCESSORY_MAX_FREQ)

    return bin_subset, mask_acc, ordered_samples

def run_pangenome_summary(meta_df: pd.DataFrame, matrix_df: pd.DataFrame, outdir: str):
    subdir = os.path.join(outdir, "Pangenome_Summary")
    make_dir(subdir)

    bin_subset, _, _ = build_binary_matrix_for_task(meta_df, matrix_df)
    if bin_subset is None:
        return pd.DataFrame(), pd.DataFrame()

    rows = []
    accum_rows = []

    rng = np.random.default_rng(RNG_SEED)

    for clade in sorted(meta_df["Lineage"].dropna().unique()):
        samples_c = meta_df.loc[meta_df["Lineage"] == clade, "Sample"].tolist()
        n = len(samples_c)
        if n == 0:
            continue

        submat = bin_subset[samples_c]
        gene_counts = submat.sum(axis=1)
        pct = gene_counts / n * 100.0

        total_seen = int((gene_counts > 0).sum())
        core_n = int((pct >= CORE_PCT).sum())
        shell_n = int(((pct >= SHELL_PCT) & (pct < CORE_PCT)).sum())
        cloud_n = int((pct < SHELL_PCT).sum())

        rows.append({
            "Clade": clade,
            "Genomes": n,
            "Total_Genes_Observed": total_seen,
            "Core_Genes": core_n,
            "Shell_Genes": shell_n,
            "Cloud_Genes": cloud_n,
            "Core_Threshold_Pct": CORE_PCT,
            "Shell_Min_Pct": SHELL_PCT
        })

        # accumulation / openness
        max_k = n
        reps = min(ACCUM_REPS, 500)
        for k in range(1, max_k + 1):
            total_genes_vals = []
            core_genes_vals = []
            for _ in range(reps):
                chosen = rng.choice(samples_c, size=k, replace=False)
                tmp = submat[list(chosen)]
                gc = tmp.sum(axis=1)
                total_genes_vals.append(int((gc > 0).sum()))
                core_genes_vals.append(int((gc == k).sum()))
            accum_rows.append({
                "Clade": clade,
                "Genomes_Sampled": k,
                "Total_Genes_Mean": float(np.mean(total_genes_vals)),
                "Total_Genes_SD": float(np.std(total_genes_vals, ddof=0)),
                "Core_Genes_Mean": float(np.mean(core_genes_vals)),
                "Core_Genes_SD": float(np.std(core_genes_vals, ddof=0))
            })

    summary_df = pd.DataFrame(rows)
    accum_df = pd.DataFrame(accum_rows)

    # openness fit
    openness_rows = []
    for clade, sub in accum_df.groupby("Clade"):
        alpha, intercept = fit_heaps_alpha(sub)
        openness_rows.append({
            "Clade": clade,
            "Heaps_Alpha": alpha,
            "Intercept": intercept
        })
    openness_df = pd.DataFrame(openness_rows)

    save_tsv_csv(summary_df, os.path.join(subdir, "Clade_Pangenome_Summary"))
    save_tsv_csv(accum_df, os.path.join(subdir, "Clade_Accumulation_Curves"))
    save_tsv_csv(openness_df, os.path.join(subdir, "Clade_Openness_Summary"))

    plot_pangenome_categories(
        summary_df,
        os.path.join(subdir, "Clade_Pangenome_StackedBar.png"),
        os.path.join(subdir, "Clade_Pangenome_StackedBar.pdf"),
        "Clade-specific pangenome composition"
    )
    plot_accumulation_curves(
        accum_df,
        os.path.join(subdir, "Clade_Accumulation_Curves.png"),
        os.path.join(subdir, "Clade_Accumulation_Curves.pdf"),
        "Clade-specific pangenome accumulation curves"
    )
    plot_core_decay_curves(
        accum_df,
        os.path.join(subdir, "Clade_Core_Decay_Curves.png"),
        os.path.join(subdir, "Clade_Core_Decay_Curves.pdf"),
        "Clade-specific core genome decay curves"
    )

    return summary_df, openness_df

def run_clade_vs_rest(meta_df: pd.DataFrame, matrix_df: pd.DataFrame, gene_info: pd.DataFrame, outdir: str):
    subdir = os.path.join(outdir, "Clade_vs_Rest")
    make_dir(subdir)

    bin_subset, mask_acc, ordered_samples = build_binary_matrix_for_task(meta_df, matrix_df)
    if bin_subset is None:
        return pd.DataFrame(), {}

    clade_sizes = meta_df["Lineage"].value_counts().to_dict()
    clades = sorted(clade_sizes.keys())
    M = len(ordered_samples)

    gene_names = matrix_df.loc[mask_acc, "Gene"].values
    annotations = matrix_df.loc[mask_acc, "Annotation"].values if "Annotation" in matrix_df.columns else np.array([""] * int(mask_acc.sum()))

    acc_np = bin_subset.loc[mask_acc].to_numpy(dtype=np.uint8)
    total_present = acc_np.sum(axis=1)

    result_rows = []
    definers_dict = {}

    for clade in clades:
        focal_samples = meta_df.loc[meta_df["Lineage"] == clade, "Sample"].tolist()
        other_samples = meta_df.loc[meta_df["Lineage"] != clade, "Sample"].tolist()

        n1 = len(focal_samples)
        n0 = len(other_samples)
        if n1 == 0 or n0 == 0:
            continue

        idx1 = [i for i, s in enumerate(ordered_samples) if s in set(focal_samples)]
        idx0 = [i for i, s in enumerate(ordered_samples) if s in set(other_samples)]

        a = acc_np[:, idx1].sum(axis=1)  # focal present
        c = acc_np[:, idx0].sum(axis=1)  # rest present
        b = n1 - a                        # focal absent
        d = n0 - c                        # rest absent

        pvals = np.array([
            fisher_exact([[int(a[i]), int(b[i])], [int(c[i]), int(d[i])]], alternative="two-sided")[1]
            for i in range(len(a))
        ])
        _, fdr, _, _ = multipletests(pvals, method="fdr_bh")

        odds = np.array([
            ((a[i] + 0.5) / (b[i] + 0.5)) / ((c[i] + 0.5) / (d[i] + 0.5))
            for i in range(len(a))
        ])

        res = pd.DataFrame({
            "Clade": clade,
            "Contrast": f"{clade} vs non-{clade}",
            "Gene": gene_names,
            "Annotation": annotations,
            "Count_Focal": a,
            "Count_Rest": c,
            "N_Focal": n1,
            "N_Rest": n0,
            "Percent_Focal": (a / n1) * 100.0,
            "Percent_Rest": (c / n0) * 100.0,
            "Odds_Ratio_0.5": odds,
            "p_value": pvals,
            "FDR": fdr,
            "Count_Total": total_present
        })
        result_rows.append(res)

        definers = res[
            (res["FDR"] <= FDR_THRESHOLD) &
            (res["Percent_Focal"] >= CLADEDEF_IN_PCT) &
            (res["Percent_Rest"] <= CLADEDEF_OUT_PCT)
        ].copy()

        definers = definers.merge(gene_info, on="Gene", how="left")
        definers_dict[clade] = definers

        if not definers.empty:
            safe = clade.replace(" ", "_")
            save_tsv_csv(definers, os.path.join(subdir, f"CladeDef_{safe}"))

    all_res = pd.concat(result_rows, ignore_index=True) if result_rows else pd.DataFrame()
    if not all_res.empty:
        all_res = all_res.merge(gene_info, on="Gene", how="left")

    summary_rows = []
    for clade in clades:
        sub = all_res[(all_res["Clade"] == clade) & (all_res["FDR"] <= FDR_THRESHOLD)] if not all_res.empty else pd.DataFrame()
        definers = definers_dict.get(clade, pd.DataFrame())
        summary_rows.append({
            "Clade": clade,
            "Samples_in_Clade": clade_sizes.get(clade, 0),
            "Accessory_Genes_Tested": int(mask_acc.sum()),
            "FDR_Significant_Genes": int(len(sub)),
            "Strict_Clade_Definers": int(len(definers)),
            "Definer_In_Pct_Cutoff": CLADEDEF_IN_PCT,
            "Definer_Out_Pct_Cutoff": CLADEDEF_OUT_PCT
        })
    summary_df = pd.DataFrame(summary_rows)

    save_tsv_csv(all_res, os.path.join(subdir, "Clade_vs_Rest_Results"))
    save_tsv_csv(summary_df, os.path.join(subdir, "Clade_vs_Rest_Summary"))

    return all_res, definers_dict

def run_pairwise(meta_df: pd.DataFrame, matrix_df: pd.DataFrame, gene_info: pd.DataFrame, outdir: str):
    subdir = os.path.join(outdir, "Pairwise")
    make_dir(subdir)

    bin_subset, mask_acc, ordered_samples = build_binary_matrix_for_task(meta_df, matrix_df)
    if bin_subset is None:
        return pd.DataFrame()

    clades = sorted(meta_df["Lineage"].dropna().unique())
    gene_names = matrix_df.loc[mask_acc, "Gene"].values
    annotations = matrix_df.loc[mask_acc, "Annotation"].values if "Annotation" in matrix_df.columns else np.array([""] * int(mask_acc.sum()))
    acc_np = bin_subset.loc[mask_acc].to_numpy(dtype=np.uint8)

    rows = []

    for i in range(len(clades)):
        for j in range(i + 1, len(clades)):
            c1 = clades[i]
            c2 = clades[j]

            s1 = meta_df.loc[meta_df["Lineage"] == c1, "Sample"].tolist()
            s2 = meta_df.loc[meta_df["Lineage"] == c2, "Sample"].tolist()
            n1, n2 = len(s1), len(s2)
            if n1 == 0 or n2 == 0:
                continue

            idx1 = [k for k, s in enumerate(ordered_samples) if s in set(s1)]
            idx2 = [k for k, s in enumerate(ordered_samples) if s in set(s2)]

            a = acc_np[:, idx1].sum(axis=1)
            c = acc_np[:, idx2].sum(axis=1)
            b = n1 - a
            d = n2 - c

            pvals = np.array([
                fisher_exact([[int(a[m]), int(b[m])], [int(c[m]), int(d[m])]], alternative="two-sided")[1]
                for m in range(len(a))
            ])
            _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
            odds = np.array([
                ((a[m] + 0.5) / (b[m] + 0.5)) / ((c[m] + 0.5) / (d[m] + 0.5))
                for m in range(len(a))
            ])

            res = pd.DataFrame({
                "Contrast": f"{c1} vs {c2}",
                "Clade1": c1,
                "Clade2": c2,
                "Gene": gene_names,
                "Annotation": annotations,
                "Count_Clade1": a,
                "Count_Clade2": c,
                "N_Clade1": n1,
                "N_Clade2": n2,
                "Percent_Clade1": (a / n1) * 100.0,
                "Percent_Clade2": (c / n2) * 100.0,
                "Odds_Ratio_0.5": odds,
                "p_value": pvals,
                "FDR": fdr
            })
            rows.append(res)

            safe = f"{c1}_vs_{c2}".replace(" ", "_")
            sig = res[res["FDR"] <= FDR_THRESHOLD].copy().merge(gene_info, on="Gene", how="left")
            save_tsv_csv(sig, os.path.join(subdir, f"Pairwise_{safe}_FDRsig"))

    all_pair = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if not all_pair.empty:
        all_pair = all_pair.merge(gene_info, on="Gene", how="left")

    save_tsv_csv(all_pair, os.path.join(subdir, "Pairwise_All_Results"))
    return all_pair

# -------------------------
# WORKER FUNCTION
# -------------------------
def run_analysis_task(args):
    label, meta_df, matrix_df, gene_info = args

    print(f"⚡ [Process {os.getpid()}] Starting: {label}")
    outdir = os.path.join(OUTROOT, label)
    os.makedirs(outdir, exist_ok=True)

    # 1. Pangenome summary + openness
    pangenome_summary_df, openness_df = run_pangenome_summary(meta_df, matrix_df, outdir)

    # 2. Clade vs rest
    clade_vs_rest_df, definers_dict = run_clade_vs_rest(meta_df, matrix_df, gene_info, outdir)

    # 3. Pairwise follow-up
    pairwise_df = run_pairwise(meta_df, matrix_df, gene_info, outdir)

    # Save reduced gene info always
    gene_info_small = gene_info.copy()
    if not clade_vs_rest_df.empty and "Gene" in clade_vs_rest_df.columns and "Gene" in gene_info.columns:
        used = set(clade_vs_rest_df["Gene"].dropna().astype(str).unique())
        gene_info_small = gene_info[gene_info["Gene"].astype(str).isin(used)].copy()
    gene_info_small = gene_info_small.drop_duplicates(subset=["Gene"]) if "Gene" in gene_info_small.columns else gene_info_small.drop_duplicates()
    save_tsv_csv(gene_info_small, os.path.join(outdir, "Gene_Info_Used"))

    # Excel
    sheet_map = {
        "Pangenome_Summary": pangenome_summary_df,
        "Openness_Summary": openness_df,
        "Clade_vs_Rest": clade_vs_rest_df,
        "Pairwise": pairwise_df,
        "Gene_Info_Used": gene_info_small
    }
    for clade_name, dfc in definers_dict.items():
        sheet_map[f"Def_{clade_name.replace(' ', '_')}"] = dfc
    save_excel_safe(outdir, sheet_map)

    return f"✅ {label}: Completed."

# -------------------------
# MAIN EXECUTION
# -------------------------
if __name__ == "__main__":
    start_time = time.time()
    print("🚀 Initializing robust Panaroo by-clade mode...")

    # 1. Data Loading
    print("📂 Loading data...")
    gene_info = load_gene_info()

    meta_full = pd.read_csv(METADATA_PATH, sep="\t", dtype=str)
    meta_full.columns = meta_full.columns.str.strip()

    if "ID" in meta_full.columns and "Sample" not in meta_full.columns:
        meta_full = meta_full.rename(columns={"ID": "Sample"})

    if "Sample" not in meta_full.columns or "Lineage" not in meta_full.columns:
        raise SystemExit("[ERROR] Clade_metadata.txt must contain Sample and Lineage")

    meta_full["Sample"] = meta_full["Sample"].map(normalize_sample_id)
    meta_full = meta_full[meta_full["Lineage"] != REFERENCE_LABEL].copy()

    required_lineages = {"Clade 1A", "Clade 1B", "Clade 2", "Clade 3", "Unassigned"}
    meta_full["Lineage"] = meta_full["Lineage"].where(
        meta_full["Lineage"].isin(required_lineages),
        other="Unassigned"
    )
    meta_full = meta_full.drop_duplicates(subset=["Sample"]).copy()

    matrix_full = pd.read_csv(PA_MATRIX_PATH, engine="pyarrow")
    if "Gene" not in matrix_full.columns:
        raise SystemExit("[ERROR] gene_presence_absence.csv missing required 'Gene' column")

    if "Annotation" not in matrix_full.columns:
        matrix_full["Annotation"] = ""

    # 2. Prepare tasks
    tasks = []

    for scheme in ["Combined_Clade1", "Split_Clade1A1B"]:
        scheme_meta = recode_lineages(meta_full, scheme)

        for mode in ["With_Unassigned", "Without_Unassigned"]:
            clade_order, include_unassigned = get_scheme_config(scheme, mode)

            mode_meta = scheme_meta.copy()
            if not include_unassigned:
                mode_meta = mode_meta[mode_meta["Lineage"] != UNASSIGNED_LABEL].copy()

            mode_meta = mode_meta[mode_meta["Lineage"].isin(clade_order)].copy()
            tasks.append((
                f"{scheme}/{mode}",
                mode_meta.copy(),
                matrix_full,
                gene_info
            ))

    print(f"🔥 Distributing {len(tasks)} scenarios across available cores...")

    # 3. Parallel execution
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(run_analysis_task, tasks))

    # 4. Final summary
    print("\n" + "\n".join(results))
    print(f"🏁 Total Execution Time: {time.time() - start_time:.2f} seconds")
