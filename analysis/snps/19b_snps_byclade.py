#!/usr/bin/env python3
import os
import re
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# =============================================================================
# CONFIGURATION
# =============================================================================
WD               = "/home/samuelajulo/SENBio/Final"
GUBBINS_DIR      = f"{WD}/gubbin_out"
METADATA_FILE    = f"{WD}/ITOL/Clade_metadata.txt"
OUTPUT_ROOT      = f"{WD}/iqtree_final/Snps_clade_integrated"
MASTER_REPORT    = f"{WD}/Snippy_output/SENBIO_RECOVERED_REPORT.csv"

# Statistics thresholds
FDR_ALPHA        = 0.05
DEF_IN_CLADE     = 40.0   # % presence in clade
DEF_OUTSIDE      = 1.0    # % allowed outside clade
UNASSIGNED_LABEL = "Unassigned"

# Expected split labels in metadata
VALID_LINEAGES_SPLIT = {"Clade 1A", "Clade 1B", "Clade 2", "Clade 3", "Unassigned"}

# Coverage bins for strict definers
# NOTE: final bin is >=90 so exact 90.0 is not lost
COVERAGE_BIN_ORDER = ["40-<60%", "60-<80%", "80-<90%", ">=90%"]

# =============================================================================
# HELPERS
# =============================================================================
def normalize_sample_id(x):
    if pd.isna(x):
        return x
    s = str(x).strip()
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

def generate_panther_list(definers_df, output_path):
    """Extract unique gene names for functional analysis."""
    if definers_df.empty:
        return 0
    genes = definers_df["GENE"].dropna().unique()
    clean_genes = [
        str(g).strip()
        for g in genes
        if str(g).lower() not in {
            "intergenic", "nan", "", "none", ".", "-", "hypothetical protein"
        }
    ]
    if clean_genes:
        with open(output_path, "w") as f:
            f.write("\n".join(sorted(set(clean_genes))))
    return len(set(clean_genes))

def safe_lower_series(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.lower()

def ensure_required_columns(df: pd.DataFrame, required_cols):
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"[ERROR] Missing required columns: {missing}")

def sanitize_clade_name(x: str) -> str:
    return str(x).replace(" ", "_").replace("/", "_")

def assign_coverage_bin(percent_in):
    if pd.isna(percent_in):
        return pd.NA
    v = float(percent_in)
    if 40.0 <= v < 60.0:
        return "40-<60%"
    elif 60.0 <= v < 80.0:
        return "60-<80%"
    elif 80.0 <= v < 90.0:
        return "80-<90%"
    elif v >= 90.0:
        return ">=90%"
    return pd.NA

def save_csv_and_xlsx(df: pd.DataFrame, csv_path: str, xlsx_path: str, sheet_name: str = "Summary"):
    df.to_csv(csv_path, index=False)
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)

def save_multi_sheet_xlsx(sheet_map: dict, xlsx_path: str):
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for sheet_name, sheet_df in sheet_map.items():
            safe_name = str(sheet_name)[:31]
            sheet_df.to_excel(writer, index=False, sheet_name=safe_name)

def build_coverage_summary_tables(summary_long: pd.DataFrame):
    if summary_long.empty:
        wide = pd.DataFrame(columns=[
            "Clade", "Samples_in_Clade", "Total_Significant_SNPs", "Total_Strict_Definers"
        ] + COVERAGE_BIN_ORDER)
        return summary_long, wide

    summary_long = summary_long.copy()
    summary_long["Coverage_Bin"] = pd.Categorical(
        summary_long["Coverage_Bin"],
        categories=COVERAGE_BIN_ORDER,
        ordered=True
    )
    summary_long = summary_long.sort_values(["Clade", "Coverage_Bin"]).reset_index(drop=True)

    wide = summary_long.pivot_table(
        index=["Clade", "Samples_in_Clade", "Total_Significant_SNPs", "Total_Strict_Definers"],
        columns="Coverage_Bin",
        values="Definer_SNP_Count",
        aggfunc="sum",
        fill_value=0
    ).reset_index()

    wide.columns.name = None

    for c in COVERAGE_BIN_ORDER:
        if c not in wide.columns:
            wide[c] = 0

    wide = wide[["Clade", "Samples_in_Clade", "Total_Significant_SNPs", "Total_Strict_Definers"] + COVERAGE_BIN_ORDER]
    return summary_long, wide

# =============================================================================
# CORE ANALYSIS
# =============================================================================
def run_stats_pipeline(meta_df, anno_df, run_outdir):
    """
    Performs hypergeometric enrichment.
    Creates:
      1. Significant_All_<clade>.csv
      2. Definer_Strict_<clade>.csv
      3. Panther_Significant_<clade>.txt
      4. Panther_Strict_<clade>.txt
      5. Detailed_Feature_Report.csv
      6. NEW: Definer coverage summaries (.csv + .xlsx)
    Returns:
      summary_long DataFrame for aggregation in main()
    """
    os.makedirs(run_outdir, exist_ok=True)

    meta_df = meta_df.copy()
    anno_df = anno_df.copy()

    meta_df["Sample"] = meta_df["Sample"].map(normalize_sample_id)
    anno_df["Sample"] = anno_df["Sample"].map(normalize_sample_id)

    meta_df = meta_df.dropna(subset=["Sample", "Lineage"]).copy()
    meta_df = meta_df.drop_duplicates(subset=["Sample"]).copy()

    subset_anno = anno_df[anno_df["Sample"].isin(meta_df["Sample"])].copy()
    if subset_anno.empty:
        print(f"[WARN] No SNPs overlap metadata for: {run_outdir}")
        pd.DataFrame().to_csv(os.path.join(run_outdir, "Detailed_Feature_Report.csv"), index=False)

        empty_long = pd.DataFrame(columns=[
            "Clade", "Samples_in_Clade", "Total_Significant_SNPs",
            "Total_Strict_Definers", "Coverage_Bin", "Definer_SNP_Count"
        ])
        save_csv_and_xlsx(
            empty_long,
            os.path.join(run_outdir, "Definer_Coverage_Summary.csv"),
            os.path.join(run_outdir, "Definer_Coverage_Summary.xlsx"),
            sheet_name="Long"
        )
        return empty_long

    # One sample contributes at most one count per SNP feature
    subset_anno = subset_anno.drop_duplicates(subset=["Sample", "POS", "REF", "ALT"]).copy()

    pivot_df = subset_anno.merge(meta_df[["Sample", "Lineage"]], on="Sample", how="inner")
    if pivot_df.empty:
        print(f"[WARN] No merged SNP+metadata rows for: {run_outdir}")
        pd.DataFrame().to_csv(os.path.join(run_outdir, "Detailed_Feature_Report.csv"), index=False)

        empty_long = pd.DataFrame(columns=[
            "Clade", "Samples_in_Clade", "Total_Significant_SNPs",
            "Total_Strict_Definers", "Coverage_Bin", "Definer_SNP_Count"
        ])
        save_csv_and_xlsx(
            empty_long,
            os.path.join(run_outdir, "Definer_Coverage_Summary.csv"),
            os.path.join(run_outdir, "Definer_Coverage_Summary.xlsx"),
            sheet_name="Long"
        )
        return empty_long

    clade_sizes = meta_df["Lineage"].value_counts().to_dict()
    clades = [c for c in meta_df["Lineage"].dropna().unique().tolist()]
    clades = sorted(clades)
    M = len(meta_df)

    counts = pivot_df.groupby(["POS", "REF", "ALT", "Lineage"]).size().unstack(fill_value=0)
    for c in clades:
        if c not in counts.columns:
            counts[c] = 0

    counts["K_total"] = counts[clades].sum(axis=1)

    anno_lookup_cols = [c for c in ["POS", "REF", "ALT", "GENE", "EFFECT", "AA_CHANGE"] if c in anno_df.columns]
    anno_lookup = anno_df[anno_lookup_cols].drop_duplicates(subset=["POS", "REF", "ALT"])

    feature_summary = []
    coverage_summary_rows = []

    for c in clades:
        n = clade_sizes[c]
        in_counts = counts[c].values
        total_k = counts["K_total"].values

        p_vals = hypergeom.sf(in_counts - 1, M, total_k, n)
        fdr = multipletests(p_vals, method="fdr_bh")[1]

        res = counts.reset_index().merge(anno_lookup, on=["POS", "REF", "ALT"], how="left")
        res["fdr"] = fdr
        res["percent_in"] = (in_counts / n) * 100.0 if n > 0 else 0.0
        res["percent_out"] = ((total_k - in_counts) / (M - n) * 100.0) if (M - n) > 0 else 0.0

        clade_clean = sanitize_clade_name(c)

        significant = res[res["fdr"] <= FDR_ALPHA].copy()
        significant.to_csv(os.path.join(run_outdir, f"Significant_All_{clade_clean}.csv"), index=False)
        p_count_sig = generate_panther_list(
            significant,
            os.path.join(run_outdir, f"Panther_Significant_{clade_clean}.txt")
        )

        definers = significant[
            (significant["percent_in"] >= DEF_IN_CLADE) &
            (significant["percent_out"] <= DEF_OUTSIDE)
        ].copy()

        definers["Coverage_Bin"] = definers["percent_in"].apply(assign_coverage_bin)

        definers.to_csv(os.path.join(run_outdir, f"Definer_Strict_{clade_clean}.csv"), index=False)
        p_count_strict = generate_panther_list(
            definers,
            os.path.join(run_outdir, f"Panther_Strict_{clade_clean}.txt")
        )

        # NEW: coverage-bin counts for strict definers
        cov_counts = definers["Coverage_Bin"].value_counts(dropna=False).to_dict()
        total_strict = len(definers)
        total_sig = len(significant)

        for bin_label in COVERAGE_BIN_ORDER:
            coverage_summary_rows.append({
                "Clade": c,
                "Samples_in_Clade": n,
                "Total_Significant_SNPs": total_sig,
                "Total_Strict_Definers": total_strict,
                "Coverage_Bin": bin_label,
                "Definer_SNP_Count": int(cov_counts.get(bin_label, 0))
            })

        intergenic_mask = safe_lower_series(definers.get("GENE", pd.Series(dtype=object))) == "intergenic"
        intergenic_count = int(intergenic_mask.sum())
        effect_counts = definers["EFFECT"].value_counts().to_dict() if "EFFECT" in definers.columns else {}

        report_entry = {
            "Clade": c,
            "Samples_in_Clade": n,
            "Total_Significant_SNPs": len(significant),
            "Strict_Defining_SNPs": len(definers),
            "Definer_Intergenic": intergenic_count,
            "Definer_Coding": len(definers) - intergenic_count,
            "Unique_Genes_Strict": p_count_strict,
            "Unique_Genes_Significant": p_count_sig
        }
        report_entry.update({f"Definer_{k}": v for k, v in effect_counts.items()})
        feature_summary.append(report_entry)

    # Existing detailed report
    pd.DataFrame(feature_summary).fillna(0).to_csv(
        os.path.join(run_outdir, "Detailed_Feature_Report.csv"),
        index=False
    )

    # NEW: per-run coverage summary outputs
    summary_long = pd.DataFrame(coverage_summary_rows)
    summary_long, summary_wide = build_coverage_summary_tables(summary_long)

    summary_csv = os.path.join(run_outdir, "Definer_Coverage_Summary.csv")
    summary_xlsx = os.path.join(run_outdir, "Definer_Coverage_Summary.xlsx")
    summary_wide_csv = os.path.join(run_outdir, "Definer_Coverage_Summary_Wide.csv")

    summary_long.to_csv(summary_csv, index=False)
    summary_wide.to_csv(summary_wide_csv, index=False)

    save_multi_sheet_xlsx(
        {
            "Long": summary_long,
            "Wide": summary_wide
        },
        summary_xlsx
    )

    return summary_long

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("🚀 Starting Integrated Clade Analysis...")

    if not os.path.exists(MASTER_REPORT):
        print(f"❌ Error: Master report not found at {MASTER_REPORT}")
        return
    if not os.path.exists(METADATA_FILE):
        print(f"❌ Error: Metadata file not found at {METADATA_FILE}")
        return

    df = pd.read_csv(MASTER_REPORT, low_memory=False)
    meta = pd.read_csv(METADATA_FILE, sep="\t", low_memory=False)

    ensure_required_columns(df, ["Sample", "POS"])
    ensure_required_columns(meta, ["Sample", "Lineage"])

    df["Sample"] = df["Sample"].map(normalize_sample_id)
    meta["Sample"] = meta["Sample"].map(normalize_sample_id)

    meta = meta.dropna(subset=["Sample", "Lineage"]).copy()
    meta = meta[meta["Lineage"] != "Reference"].copy()
    meta["Lineage"] = meta["Lineage"].where(meta["Lineage"].isin(VALID_LINEAGES_SPLIT), other="Unassigned")
    meta = meta.drop_duplicates(subset=["Sample"]).copy()

    # Ensure core columns exist for downstream export
    for col in ["REF", "ALT", "GENE", "EFFECT", "AA_CHANGE"]:
        if col not in df.columns:
            df[col] = pd.NA

    # Numeric POS
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df = df.dropna(subset=["Sample", "POS"]).copy()
    df["POS"] = df["POS"].astype(int)

    # Recombination filtering
    df["is_recombination"] = False
    gff_path = os.path.join(GUBBINS_DIR, "senbio_clean.recombination_predictions.gff")

    if os.path.exists(gff_path):
        print(f"🧬 Flagging recombination from GFF: {gff_path}")
        with open(gff_path, "r") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9:
                    continue
                start, end = int(p[3]), int(p[4])
                attributes = p[8]
                if 'taxa="' in attributes:
                    taxa_str = attributes.split('taxa="')[1].split('"')[0]
                    taxa_list = taxa_str.replace(";", " ").strip().split()
                    taxa_list = [normalize_sample_id(x) for x in taxa_list]
                    df.loc[
                        (df["Sample"].isin(taxa_list)) &
                        (df["POS"] >= start) &
                        (df["POS"] <= end),
                        "is_recombination"
                    ] = True
        print(f"✅ Flagged {int(df['is_recombination'].sum())} SNP records as recombination.")
    else:
        print("⚠️ GFF NOT FOUND. Proceeding with 0 recombination flags.")

    scenarios = [
        ("all_snps", df.copy()),
        ("vertical_only", df[df["is_recombination"] == False].copy())
    ]

    schemes = ["Combined_Clade1", "Split_Clade1A1B"]
    modes = ["with_unassigned", "without_unassigned"]

    all_run_summaries = []

    for scenario_label, current_data in scenarios:
        for scheme in schemes:
            scheme_meta = recode_lineages(meta, scheme)

            for mode in modes:
                if mode == "with_unassigned":
                    filtered_meta = scheme_meta.copy()
                else:
                    filtered_meta = scheme_meta[scheme_meta["Lineage"] != UNASSIGNED_LABEL].copy()

                run_outdir = os.path.join(
                    OUTPUT_ROOT,
                    scheme,
                    scenario_label,
                    mode
                )

                print(f"🔵 Running Scenario={scenario_label} | Scheme={scheme} | Mode={mode}")
                run_summary = run_stats_pipeline(filtered_meta, current_data, run_outdir)

                if run_summary is not None and not run_summary.empty:
                    run_summary = run_summary.copy()
                    run_summary["Scheme"] = scheme
                    run_summary["Scenario"] = scenario_label
                    run_summary["Mode"] = mode
                    all_run_summaries.append(run_summary)

    # NEW: overall combined summaries across all run outputs
    if all_run_summaries:
        overall_long = pd.concat(all_run_summaries, ignore_index=True)
        overall_long["Coverage_Bin"] = pd.Categorical(
            overall_long["Coverage_Bin"],
            categories=COVERAGE_BIN_ORDER,
            ordered=True
        )
        overall_long = overall_long.sort_values(
            ["Scheme", "Scenario", "Mode", "Clade", "Coverage_Bin"]
        ).reset_index(drop=True)

        overall_csv = os.path.join(OUTPUT_ROOT, "Overall_Definer_Coverage_Summary.csv")
        overall_xlsx = os.path.join(OUTPUT_ROOT, "Overall_Definer_Coverage_Summary.xlsx")
        overall_wide_csv = os.path.join(OUTPUT_ROOT, "Overall_Definer_Coverage_Summary_Wide.csv")
        overall_by_clade_csv = os.path.join(OUTPUT_ROOT, "Overall_Definer_Coverage_ByClade.csv")

        overall_wide = overall_long.pivot_table(
            index=["Scheme", "Scenario", "Mode", "Clade", "Samples_in_Clade", "Total_Significant_SNPs", "Total_Strict_Definers"],
            columns="Coverage_Bin",
            values="Definer_SNP_Count",
            aggfunc="sum",
            fill_value=0
        ).reset_index()
        overall_wide.columns.name = None

        for c in COVERAGE_BIN_ORDER:
            if c not in overall_wide.columns:
                overall_wide[c] = 0

        # overall summary for each clade across all created outputs
        overall_by_clade = overall_long.groupby(
            ["Clade", "Coverage_Bin"],
            as_index=False
        ).agg(
            Definer_SNP_Count=("Definer_SNP_Count", "sum"),
            Runs_Contributing=("Scheme", "count")
        )

        overall_by_clade["Coverage_Bin"] = pd.Categorical(
            overall_by_clade["Coverage_Bin"],
            categories=COVERAGE_BIN_ORDER,
            ordered=True
        )
        overall_by_clade = overall_by_clade.sort_values(["Clade", "Coverage_Bin"]).reset_index(drop=True)

        overall_long.to_csv(overall_csv, index=False)
        overall_wide.to_csv(overall_wide_csv, index=False)
        overall_by_clade.to_csv(overall_by_clade_csv, index=False)

        overall_by_clade_wide = overall_by_clade.pivot_table(
            index="Clade",
            columns="Coverage_Bin",
            values="Definer_SNP_Count",
            aggfunc="sum",
            fill_value=0
        ).reset_index()
        overall_by_clade_wide.columns.name = None

        for c in COVERAGE_BIN_ORDER:
            if c not in overall_by_clade_wide.columns:
                overall_by_clade_wide[c] = 0

        save_multi_sheet_xlsx(
            {
                "By_Run_Long": overall_long,
                "By_Run_Wide": overall_wide,
                "By_Clade_Long": overall_by_clade,
                "By_Clade_Wide": overall_by_clade_wide
            },
            overall_xlsx
        )

    print(f"\n🎉 ANALYSIS COMPLETE. Results located in: {OUTPUT_ROOT}")

if __name__ == "__main__":
    main()
