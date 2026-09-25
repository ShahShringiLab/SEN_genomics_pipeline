#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import itertools
import warnings

# --- 1) DEPENDENCY MANAGER ---
def check_and_install_dependencies():
    required = {
        "pandas": "pandas",
        "numpy": "numpy",
        "seaborn": "seaborn",
        "matplotlib": "matplotlib",
        "scipy": "scipy",
        "statsmodels": "statsmodels",
        "openpyxl": "openpyxl",
        "PIL": "Pillow",
    }
    for imp, pip_name in required.items():
        try:
            __import__(imp)
        except ImportError:
            subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name])

check_and_install_dependencies()

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
# CONFIG
# =============================================================================
BASE_DIR = "/home/samuelajulo/SENBio/Final"

vfdb_report_path = os.path.join(BASE_DIR, "abricate_results", "vfdb", "master_vfdb_report.tsv")

metadata_path_candidates = [
    os.path.join(BASE_DIR, "iqtree_final", "ITOL", "Clade_metadata.txt"),
    os.path.join(BASE_DIR, "ITOL", "Clade_metadata.txt"),
    os.path.join(BASE_DIR, "iqtree", "ITOL", "Clade_metadata.txt"),
]
metadata_path = next((p for p in metadata_path_candidates if os.path.exists(p)), None)
if metadata_path is None:
    raise SystemExit("[ERROR] Could not find Clade_metadata.txt in expected locations.")

ROOT_OUT_DIR = os.path.join(BASE_DIR, "iqtree_final", "VFDB_clade")
os.makedirs(ROOT_OUT_DIR, exist_ok=True)

# Presence definition
MIN_COVERAGE = 80.0
MIN_IDENTITY = 90.0

# Stats
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

# Heatmap colors (0/1)
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
    s = s.strip("_")
    s = re.sub(r"(_trimmed.*)$", "", s)
    s = re.sub(r"(_pure.*)$", "", s)
    s = re.sub(r"(_R[12])$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"(_contigs)$", "", s, flags=re.IGNORECASE)
    m = re.search(r"(SRR\d+)", s)
    return m.group(1) if m else s

def extract_category_from_product(prod: str) -> str:
    """
    PRODUCT often contains two bracket blocks:
      ... [CATEGORY STUFF] [Organism/context]
    We take the FIRST [...] only and drop brackets.
    If none found, return 'Uncategorized'.
    """
    if pd.isna(prod):
        return "Uncategorized"
    s = str(prod)
    m = re.search(r"\[([^\]]+)\]", s)
    if m:
        return m.group(1).strip()
    return "Uncategorized"

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

def read_tsv_robust(path: str, badlines_out: str):
    expected_nf, bad = tsv_integrity_check(path)
    print(f"[INFO] master header fields = {expected_nf}")
    if bad:
        print(f"[WARN] Found {len(bad)} malformed lines in master TSV.")
        pd.DataFrame(bad).to_csv(badlines_out, sep="\t", index=False)
        print(f"[WARN] Bad line report saved: {badlines_out}")
        df = pd.read_csv(path, sep="\t", engine="python", on_bad_lines="skip", low_memory=False)
        print(f"[WARN] Loaded with skipped bad lines. Rows now={len(df)}")
        return df
    return pd.read_csv(path, sep="\t", low_memory=False)

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
        "excel": os.path.join(base_out, "VFDB_ByClade_Summary.xlsx"),
        "summary": os.path.join(base_out, "VFDB_ByClade_Summary.tsv"),
        "omnibus": os.path.join(base_out, "VFDB_Omnibus_Stats.tsv"),
        "pairwise": os.path.join(base_out, "VFDB_Pairwise_Fisher.tsv"),
        "catalog": os.path.join(base_out, "VFDB_Gene_Catalog.tsv"),
        "matrix": os.path.join(base_out, "VFDB_PresenceAbsence_Matrix.tsv"),
        "badlines": os.path.join(base_out, "vfdb_master_bad_lines.tsv"),
        "missing": os.path.join(base_out, "samples_in_metadata_missing_in_vfdb.tsv"),
        "heat_png": os.path.join(base_out, "VFDB_Heatmap_ByClade.png"),
        "heat_pdf": os.path.join(base_out, "VFDB_Heatmap_ByClade.pdf"),
        "heat_tif": os.path.join(base_out, "VFDB_Heatmap_ByClade.tif"),
    }

# =============================================================================
# LOAD MASTER DATA ONCE
# =============================================================================
print("📂 Loading VFDB data...")

base_out_for_badlines = os.path.join(ROOT_OUT_DIR, "_tmp")
os.makedirs(base_out_for_badlines, exist_ok=True)
tmp_badlines = os.path.join(base_out_for_badlines, "vfdb_master_bad_lines.tsv")

if not os.path.exists(vfdb_report_path):
    raise SystemExit(f"[ERROR] Missing: {vfdb_report_path}")
if not os.path.exists(metadata_path):
    raise SystemExit(f"[ERROR] Missing: {metadata_path}")

vf_df = read_tsv_robust(vfdb_report_path, tmp_badlines)
meta_df_master = pd.read_csv(metadata_path, sep="\t", low_memory=False)

if "Sample" not in meta_df_master.columns or "Lineage" not in meta_df_master.columns:
    raise SystemExit("[ERROR] Clade_metadata.txt must have columns: Sample, Lineage")

vf_df.rename(columns={vf_df.columns[0]: "Sample"} if vf_df.columns[0] != "Sample" else {}, inplace=True)

vf_df["Sample_raw"] = vf_df["Sample"].astype(str)
vf_df = vf_df[~vf_df["Sample_raw"].str.startswith(EXCLUDE_SMOKE_PREFIX)].copy()

vf_df["Sample"] = vf_df["Sample_raw"].map(normalize_sample_id)
meta_df_master["Sample"] = meta_df_master["Sample"].map(normalize_sample_id)

meta_df_master = meta_df_master[meta_df_master["Lineage"] != "Reference"].copy()

required_lineages = {"Clade 1A", "Clade 1B", "Clade 2", "Clade 3", "Unassigned"}
meta_df_master["Lineage"] = meta_df_master["Lineage"].where(
    meta_df_master["Lineage"].isin(required_lineages),
    other="Unassigned"
)

# robust column detection
colmap = {str(c).strip().lower(): c for c in vf_df.columns}

gene_col = None
for cand in ["gene"]:
    if cand in colmap:
        gene_col = colmap[cand]
        break
if gene_col is None:
    raise SystemExit("[ERROR] Could not find GENE column in master_vfdb_report.tsv")

coverage_col = None
for cand in [
    "%coverage",
    "% coverage",
    "% coverage of reference",
    "coverage"
]:
    if cand in colmap:
        coverage_col = colmap[cand]
        break

identity_col = None
for cand in [
    "%identity",
    "% identity",
    "% identity to reference",
    "identity"
]:
    if cand in colmap:
        identity_col = colmap[cand]
        break

before = len(vf_df)

vf_df = vf_df[vf_df[gene_col].notna()].copy()

if coverage_col is not None:
    vf_df[coverage_col] = pd.to_numeric(vf_df[coverage_col], errors="coerce")
    vf_df = vf_df[vf_df[coverage_col] >= MIN_COVERAGE].copy()
else:
    print("[WARN] No coverage column found in master_vfdb_report.tsv; skipping coverage filter.")

if identity_col is not None:
    vf_df[identity_col] = pd.to_numeric(vf_df[identity_col], errors="coerce")
    vf_df = vf_df[vf_df[identity_col] >= MIN_IDENTITY].copy()
else:
    print("[WARN] No identity column found in master_vfdb_report.tsv; skipping identity filter.")

after = len(vf_df)
print(f"[INFO] Applied filters: coverage >= {MIN_COVERAGE} and identity >= {MIN_IDENTITY}: {before} -> {after}")
print(f"[INFO] Coverage column used: {coverage_col}")
print(f"[INFO] Identity column used: {identity_col}")
print(f"[INFO] Gene column used: {gene_col}")
print(f"[INFO] Raw metadata lineages found: {sorted(meta_df_master['Lineage'].dropna().unique())}")

vf_df = vf_df.rename(columns={gene_col: "GENE"})

if "PRODUCT" in vf_df.columns:
    vf_df["Category"] = vf_df["PRODUCT"].map(extract_category_from_product)
else:
    vf_df["PRODUCT"] = np.nan
    vf_df["Category"] = "Uncategorized"

# =============================================================================
# RUN BOTH SCHEMES × MODES
# =============================================================================
for scheme in ["Combined_Clade1", "Split_Clade1A1B"]:
    for mode in ["With_Unassigned", "Without_Unassigned"]:
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

        vf_curr = vf_df[vf_df["Sample"].isin(set(allowed_samples))].copy()
        print(f"[INFO] VFDB rows after filters={len(vf_curr)}")

        cat_cols = ["GENE", "Category", "PRODUCT"]
        if "ACCESSION" in vf_curr.columns:
            cat_cols.append("ACCESSION")

        catalog = (
            vf_curr.groupby(cat_cols, dropna=False).size().reset_index(name="n")
            .sort_values(["GENE", "n"], ascending=[True, False])
            .drop_duplicates(subset=["GENE"], keep="first")
            .drop(columns=["n"])
        )
        catalog.to_csv(out["catalog"], sep="\t", index=False)

        print("🧱 Building presence/absence matrix...")
        pa = vf_curr.pivot_table(
            index="Sample",
            columns="GENE",
            aggfunc="size",
            fill_value=0
        )
        pa = (pa > 0).astype(int)

        missing_samples = meta_df.loc[~meta_df["Sample"].isin(pa.index), "Sample"].tolist()
        if missing_samples:
            pd.DataFrame({"Sample": missing_samples}).to_csv(out["missing"], sep="\t", index=False)
            print(f"[WARN] {len(missing_samples)} metadata samples missing in VFDB hits; they will be added as all-zero rows.")

        pa = pa.reindex(meta_df["Sample"], fill_value=0)
        pa.to_csv(out["matrix"], sep="\t")
        print(f"[OK] Matrix: {out['matrix']} shape={pa.shape}")

        print("📊 Running omnibus stats per gene...")

        clade_sizes = meta_df["Lineage"].value_counts().to_dict()

        def clade_header(c):
            return f"{c} (n={clade_sizes.get(c,0)})"

        counts_by_clade = {}
        for c in clade_order:
            samples_c = meta_df.loc[meta_df["Lineage"] == c, "Sample"].values
            counts_by_clade[c] = pa.loc[samples_c].sum(axis=0) if len(samples_c) else pd.Series(0, index=pa.columns)

        total_present = pa.sum(axis=0)
        gene_to_cat = dict(zip(catalog["GENE"], catalog["Category"])) if not catalog.empty else {}

        omnibus_rows = []
        used_clades_by_gene = {}

        for gene in pa.columns:
            pres = {c: int(counts_by_clade[c].loc[gene]) for c in clade_order}
            ns   = {c: int(clade_sizes.get(c, 0)) for c in clade_order}
            absn = {c: ns[c] - pres[c] for c in clade_order}
            tot_pres = int(total_present.loc[gene])

            variable_clades = [c for c in clade_order if ns[c] > 0 and (0 < pres[c] < ns[c])]
            used_clades_by_gene[gene] = variable_clades

            pval = np.nan
            method = ""
            chi2 = np.nan
            dof = np.nan
            reason = ""

            if tot_pres < MIN_TOTAL_PRESENT:
                reason = f"rare (<{MIN_TOTAL_PRESENT} total present)"
            elif len(variable_clades) < 2:
                reason = "insufficient variable clades after zero/full exclusion"
            else:
                table = np.array([[pres[c], absn[c]] for c in variable_clades], dtype=int)
                try:
                    pval, method, chi2, dof = chi2_with_optional_montecarlo(table)
                except Exception as e:
                    reason = f"chi2 failed: {e}"

            omnibus_rows.append({
                "GENE": gene,
                "Category": gene_to_cat.get(gene, "Uncategorized"),
                "Total_Present": tot_pres,
                "Tested_Clades": ",".join(variable_clades),
                "Test_Method": method,
                "Chi2": chi2,
                "DoF": dof,
                "Stats p-value": pval,
                "Stats FDR": np.nan,
                "Reason_Not_Tested": reason,
                **{f"{c}_present": pres[c] for c in clade_order},
                **{f"{c}_n": ns[c] for c in clade_order},
            })

        omnibus_df = pd.DataFrame(omnibus_rows)

        p_num = pd.to_numeric(omnibus_df["Stats p-value"], errors="coerce")
        tested_mask = p_num.notna()

        if tested_mask.sum() > 0:
            _, fdr_vals, _, _ = multipletests(p_num.loc[tested_mask].values.astype(float), method="fdr_bh")
            omnibus_df.loc[tested_mask, "Stats FDR"] = fdr_vals

        omnibus_df.to_csv(out["omnibus"], sep="\t", index=False)
        print(f"[OK] Omnibus stats: {out['omnibus']}")

        print("📊 Running post-hoc pairwise Fisher (Holm within gene) for significant genes...")

        fdr_num = pd.to_numeric(omnibus_df["Stats FDR"], errors="coerce")
        sig_genes = omnibus_df.loc[tested_mask & (fdr_num <= ALPHA), "GENE"].tolist()

        pairwise_rows = []
        cld_letters_by_gene = {g: {c: "" for c in clade_order} for g in pa.columns}

        for gene in sig_genes:
            pres = {c: int(counts_by_clade[c].loc[gene]) for c in clade_order}
            ns   = {c: int(clade_sizes.get(c, 0)) for c in clade_order}
            absn = {c: ns[c] - pres[c] for c in clade_order}

            variable_clades = used_clades_by_gene.get(gene, [])
            if len(variable_clades) < 2:
                continue

            pairs = list(itertools.combinations(variable_clades, 2))
            pvals = []
            tmp_rows = []

            for a, b in pairs:
                a1, a0 = pres[a], absn[a]
                b1, b0 = pres[b], absn[b]

                _, p_raw = fisher_exact([[a1, a0], [b1, b0]], alternative="two-sided")
                or05 = odds_ratio_with_haldene(a1, a0, b1, b0, add=0.5)

                tmp_rows.append({
                    "GENE": gene,
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

            cld_var = make_cld_letters(variable_clades, p_adj_map, alpha=ALPHA)
            for c in variable_clades:
                cld_letters_by_gene[gene][c] = cld_var.get(c, "")

        pairwise_df = pd.DataFrame(pairwise_rows)
        pairwise_df.to_csv(out["pairwise"], sep="\t", index=False)
        print(f"[OK] Pairwise stats: {out['pairwise']}")

        print("🧾 Building final summary table...")

        p_map = dict(zip(omnibus_df["GENE"], omnibus_df["Stats p-value"]))
        f_map = dict(zip(omnibus_df["GENE"], omnibus_df["Stats FDR"]))
        reason_map = dict(zip(omnibus_df["GENE"], omnibus_df["Reason_Not_Tested"]))
        total_map = dict(zip(omnibus_df["GENE"], omnibus_df["Total_Present"]))

        gene_to_cat = dict(zip(catalog["GENE"], catalog["Category"])) if not catalog.empty else {}
        gene_to_prod = dict(zip(catalog["GENE"], catalog["PRODUCT"])) if ("PRODUCT" in catalog.columns and not catalog.empty) else {}

        summary_rows = []
        for gene in pa.columns:
            row = {
                "Category": gene_to_cat.get(gene, "Uncategorized"),
                "GENE": gene,
                "PRODUCT": gene_to_prod.get(gene, ""),
            }
            letters = cld_letters_by_gene.get(gene, {c: "" for c in clade_order})

            for c in clade_order:
                n = int(clade_sizes.get(c, 0))
                cnt = int(counts_by_clade[c].loc[gene]) if n > 0 else 0
                row[clade_header(c)] = format_cell(cnt, n, letters.get(c, ""))

            pval = p_map.get(gene, np.nan)
            fdr  = f_map.get(gene, np.nan)
            row["Stats p-value"] = "" if pd.isna(pval) else float(pval)
            row["Stats FDR"]     = "" if pd.isna(fdr) else float(fdr)

            row["Total (n)"] = int(total_map.get(gene, int(pa[gene].sum())))
            row["Tested"] = "" if (reason_map.get(gene, "") == "") else "No"
            row["Reason_Not_Tested"] = reason_map.get(gene, "")

            summary_rows.append(row)

        summary_df = pd.DataFrame(summary_rows)
        summary_df["_cat"] = summary_df["Category"].fillna("Uncategorized").astype(str)
        summary_df["_tested_rank"] = np.where(summary_df["Reason_Not_Tested"].fillna("") == "", 0, 1)
        summary_df["_total"] = pd.to_numeric(summary_df["Total (n)"], errors="coerce").fillna(0).astype(int)

        summary_df = summary_df.sort_values(
            by=["_cat", "_total", "_tested_rank", "_total", "GENE"],
            ascending=[True, False, True, False, True]
        ).drop(columns=["_cat", "_tested_rank", "_total"]).reset_index(drop=True)

        summary_df.to_csv(out["summary"], sep="\t", index=False)
        print(f"[OK] Summary table: {out['summary']}")

        print("🎨 Rendering heatmap...")

        col_order = pa.sum(axis=0).sort_values(ascending=False).index.tolist()
        plot_data = pa[col_order].copy()

        row_colors = meta_df.set_index("Sample")["Lineage"].map(clade_palette)

        g = sns.clustermap(
            plot_data,
            cmap=my_cmap,
            row_cluster=False,
            col_cluster=True,
            row_colors=row_colors,
            figsize=(35, 90),
            yticklabels=True,
            xticklabels=True,
            linewidths=0,
            cbar_pos=None
        )

        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize=5)
        plt.setp(g.ax_heatmap.get_yticklabels(), fontsize=2)

        legend_elements = [
            Patch(facecolor="#FF0000", label="Present (1)"),
            Patch(facecolor="#FFFFFF", edgecolor="gray", label="Absent (0)"),
        ] + [Patch(facecolor=clade_palette[c], label=c) for c in clade_order if clade_sizes.get(c, 0) > 0]

        plt.legend(handles=legend_elements, title="Legend", loc="upper left", bbox_to_anchor=(1.05, 1))

        g.savefig(out["heat_pdf"], bbox_inches="tight")
        g.savefig(out["heat_png"], dpi=600, bbox_inches="tight")
        g.savefig(out["heat_tif"], dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
        plt.close("all")

        print(f"[OK] Heatmap PDF: {out['heat_pdf']}")
        print(f"[OK] Heatmap PNG: {out['heat_png']}")
        print(f"[OK] Heatmap TIF: {out['heat_tif']}")

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
            catalog.to_excel(xw, sheet_name="Gene_Catalog", index=False)
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

        if os.path.exists(out["badlines"]):
            print(f"[NOTE] vfdb master TSV had malformed lines -> {out['badlines']}")

        print("\n[INFO] Quick clade counts:")
        print(meta_df["Lineage"].value_counts().reindex(clade_order).fillna(0).astype(int))
        print(f"[INFO] Significant genes (FDR <= {ALPHA}): {len(sig_genes)}")

print("\n🎉 ALL DONE")
print(f"[OK] Results root: {ROOT_OUT_DIR}")
