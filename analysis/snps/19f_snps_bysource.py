#!/usr/bin/env python3
import os
import re
from pathlib import Path
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

# =============================================================================
# CONFIGURATION
# =============================================================================
REPO_ROOT = Path(__file__).resolve().parents[2]
WD = str(Path(os.environ.get("SEN_ROOT", REPO_ROOT)))

GUBBINS_DIR   = f"{WD}/gubbin_out"
MASTER_REPORT = f"{WD}/Snippy_output/SENBIO_RECOVERED_REPORT.csv"
CLADE_META    = os.environ.get("SEN_CLADE_METADATA", f"{WD}/metadata/final_clade_metadata.tsv")
SOURCE_FILE   = f"{WD}/itol_1_source.txt"

OUTPUT_ROOT   = f"{WD}/iqtree_final/SNPs_by_source_modest_variation"

# analysis scope
FOCAL_CLADES  = {"Clade 2", "Clade 3"}
FOCAL_SOURCES = ["Human", "Chicken-US"]

# statistics / defining thresholds
FDR_ALPHA     = 0.05
DEF_IN_GROUP  = 40.0
DEF_OUTSIDE   = 1.0

# coverage bins
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

def ensure_required_columns(df: pd.DataFrame, required_cols):
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise SystemExit(f"[ERROR] Missing required columns: {missing}")

def sanitize_name(x: str) -> str:
    return str(x).replace(" ", "_").replace("/", "_").replace("-", "_")

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

def safe_lower_series(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.lower()

def generate_panther_list(definers_df, output_path):
    if definers_df.empty:
        return 0
    genes = definers_df["GENE"].dropna().unique() if "GENE" in definers_df.columns else []
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

def save_multi_sheet_xlsx(sheet_map: dict, xlsx_path: str):
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        for sheet_name, sheet_df in sheet_map.items():
            safe_name = str(sheet_name)[:31]
            sheet_df.to_excel(writer, index=False, sheet_name=safe_name)

def build_coverage_summary_tables(summary_long: pd.DataFrame):
    if summary_long.empty:
        wide = pd.DataFrame(columns=[
            "Clade", "Source_Target", "Samples_in_Clade", "Target_Samples",
            "Total_Significant_SNPs", "Total_Strict_Definers"
        ] + COVERAGE_BIN_ORDER)
        return summary_long, wide

    summary_long = summary_long.copy()
    summary_long["Coverage_Bin"] = pd.Categorical(
        summary_long["Coverage_Bin"],
        categories=COVERAGE_BIN_ORDER,
        ordered=True
    )
    summary_long = summary_long.sort_values(
        ["Clade", "Source_Target", "Coverage_Bin"]
    ).reset_index(drop=True)

    wide = summary_long.pivot_table(
        index=[
            "Clade", "Source_Target", "Samples_in_Clade", "Target_Samples",
            "Total_Significant_SNPs", "Total_Strict_Definers"
        ],
        columns="Coverage_Bin",
        values="Definer_SNP_Count",
        aggfunc="sum",
        fill_value=0
    ).reset_index()

    wide.columns.name = None
    for c in COVERAGE_BIN_ORDER:
        if c not in wide.columns:
            wide[c] = 0

    wide = wide[[
        "Clade", "Source_Target", "Samples_in_Clade", "Target_Samples",
        "Total_Significant_SNPs", "Total_Strict_Definers"
    ] + COVERAGE_BIN_ORDER]

    return summary_long, wide

def parse_source_colorstrip(path):
    """
    Parse iTOL DATASET_COLORSTRIP file:
    sample,color,label
    """
    rows = []
    in_data = False

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("DATA"):
                in_data = True
                continue
            if not in_data:
                continue

            parts = [x.strip() for x in line.split(",")]
            if len(parts) < 3:
                continue

            sample = normalize_sample_id(parts[0])
            source = parts[2]
            rows.append((sample, source))

    df = pd.DataFrame(rows, columns=["Sample", "Source"])
    df = df.dropna(subset=["Sample", "Source"]).drop_duplicates(subset=["Sample"])
    return df

# =============================================================================
# CORE ANALYSIS
# =============================================================================
def run_source_pipeline(meta_df, anno_df, run_outdir):
    """
    Within-clade source-associated SNP analysis:
      For each focal clade (2, 3)
      For each focal source (Human, Chicken-US)
      compare target source vs rest of sources inside that clade
    """
    os.makedirs(run_outdir, exist_ok=True)

    meta_df = meta_df.copy()
    anno_df = anno_df.copy()

    meta_df["Sample"] = meta_df["Sample"].map(normalize_sample_id)
    anno_df["Sample"] = anno_df["Sample"].map(normalize_sample_id)

    meta_df = meta_df.dropna(subset=["Sample", "Lineage", "Source"]).copy()
    meta_df = meta_df.drop_duplicates(subset=["Sample"]).copy()

    subset_anno = anno_df[anno_df["Sample"].isin(meta_df["Sample"])].copy()
    if subset_anno.empty:
        print(f"[WARN] No SNPs overlap metadata for: {run_outdir}")
        pd.DataFrame().to_csv(os.path.join(run_outdir, "Detailed_Source_Report.csv"), index=False)

        empty_long = pd.DataFrame(columns=[
            "Clade", "Source_Target", "Samples_in_Clade", "Target_Samples",
            "Total_Significant_SNPs", "Total_Strict_Definers",
            "Coverage_Bin", "Definer_SNP_Count"
        ])
        empty_long.to_csv(os.path.join(run_outdir, "Source_Definer_Coverage_Summary.csv"), index=False)
        return empty_long

    subset_anno = subset_anno.drop_duplicates(subset=["Sample", "POS", "REF", "ALT"]).copy()
    pivot_df = subset_anno.merge(meta_df[["Sample", "Lineage", "Source"]], on="Sample", how="inner")

    if pivot_df.empty:
        print(f"[WARN] No merged SNP+metadata rows for: {run_outdir}")
        pd.DataFrame().to_csv(os.path.join(run_outdir, "Detailed_Source_Report.csv"), index=False)

        empty_long = pd.DataFrame(columns=[
            "Clade", "Source_Target", "Samples_in_Clade", "Target_Samples",
            "Total_Significant_SNPs", "Total_Strict_Definers",
            "Coverage_Bin", "Definer_SNP_Count"
        ])
        empty_long.to_csv(os.path.join(run_outdir, "Source_Definer_Coverage_Summary.csv"), index=False)
        return empty_long

    anno_lookup_cols = [c for c in ["POS", "REF", "ALT", "GENE", "EFFECT", "AA_CHANGE"] if c in anno_df.columns]
    anno_lookup = anno_df[anno_lookup_cols].drop_duplicates(subset=["POS", "REF", "ALT"])

    feature_summary = []
    coverage_summary_rows = []

    for clade in sorted(FOCAL_CLADES):
        clade_meta = meta_df[meta_df["Lineage"] == clade].copy()
        if clade_meta.empty:
            continue

        clade_n = len(clade_meta)

        clade_anno = pivot_df[pivot_df["Lineage"] == clade].copy()
        if clade_anno.empty:
            continue

        sources_present = set(clade_meta["Source"].unique())

        for source_target in FOCAL_SOURCES:
            if source_target not in sources_present:
                print(f"[INFO] Skipping {clade} | {source_target}: no target samples present")
                continue

            target_samples = set(clade_meta.loc[clade_meta["Source"] == source_target, "Sample"])
            rest_samples = set(clade_meta.loc[clade_meta["Source"] != source_target, "Sample"])

            n_target = len(target_samples)
            n_rest = len(rest_samples)

            if n_target == 0 or n_rest == 0:
                print(f"[INFO] Skipping {clade} | {source_target}: target={n_target}, rest={n_rest}")
                continue

            counts = clade_anno.groupby(["POS", "REF", "ALT", "Source"]).size().unstack(fill_value=0)

            # target count
            if source_target not in counts.columns:
                counts[source_target] = 0

            target_count = counts[source_target].copy()
            total_count = counts.sum(axis=1)

            p_vals = hypergeom.sf(target_count.values - 1, clade_n, total_count.values, n_target)
            fdr = multipletests(p_vals, method="fdr_bh")[1]

            res = counts.reset_index().merge(anno_lookup, on=["POS", "REF", "ALT"], how="left")
            res["fdr"] = fdr
            res["target_source"] = source_target
            res["clade"] = clade
            res["count_in_target"] = target_count.values
            res["count_total"] = total_count.values
            res["percent_in"] = (target_count.values / n_target) * 100.0 if n_target > 0 else 0.0
            res["percent_out"] = ((total_count.values - target_count.values) / n_rest * 100.0) if n_rest > 0 else 0.0

            clade_clean = sanitize_name(clade)
            source_clean = sanitize_name(source_target)

            significant = res[res["fdr"] <= FDR_ALPHA].copy()
            significant.to_csv(
                os.path.join(run_outdir, f"Significant_All_{clade_clean}_{source_clean}_vs_rest.csv"),
                index=False
            )
            p_count_sig = generate_panther_list(
                significant,
                os.path.join(run_outdir, f"Panther_Significant_{clade_clean}_{source_clean}_vs_rest.txt")
            )

            definers = significant[
                (significant["percent_in"] >= DEF_IN_GROUP) &
                (significant["percent_out"] <= DEF_OUTSIDE)
            ].copy()

            definers["Coverage_Bin"] = definers["percent_in"].apply(assign_coverage_bin)
            definers.to_csv(
                os.path.join(run_outdir, f"Definer_Strict_{clade_clean}_{source_clean}_vs_rest.csv"),
                index=False
            )
            p_count_strict = generate_panther_list(
                definers,
                os.path.join(run_outdir, f"Panther_Strict_{clade_clean}_{source_clean}_vs_rest.txt")
            )

            cov_counts = definers["Coverage_Bin"].value_counts(dropna=False).to_dict()
            total_strict = len(definers)
            total_sig = len(significant)

            for bin_label in COVERAGE_BIN_ORDER:
                coverage_summary_rows.append({
                    "Clade": clade,
                    "Source_Target": source_target,
                    "Samples_in_Clade": clade_n,
                    "Target_Samples": n_target,
                    "Total_Significant_SNPs": total_sig,
                    "Total_Strict_Definers": total_strict,
                    "Coverage_Bin": bin_label,
                    "Definer_SNP_Count": int(cov_counts.get(bin_label, 0))
                })

            intergenic_mask = safe_lower_series(definers.get("GENE", pd.Series(dtype=object))) == "intergenic"
            intergenic_count = int(intergenic_mask.sum())
            effect_counts = definers["EFFECT"].value_counts().to_dict() if "EFFECT" in definers.columns else {}

            report_entry = {
                "Clade": clade,
                "Source_Target": source_target,
                "Samples_in_Clade": clade_n,
                "Target_Samples": n_target,
                "Rest_Samples": n_rest,
                "Total_Significant_SNPs": total_sig,
                "Strict_Defining_SNPs": total_strict,
                "Definer_Intergenic": intergenic_count,
                "Definer_Coding": total_strict - intergenic_count,
                "Unique_Genes_Strict": p_count_strict,
                "Unique_Genes_Significant": p_count_sig
            }
            report_entry.update({f"Definer_{k}": v for k, v in effect_counts.items()})
            feature_summary.append(report_entry)

    detailed_df = pd.DataFrame(feature_summary).fillna(0)
    detailed_df.to_csv(os.path.join(run_outdir, "Detailed_Source_Report.csv"), index=False)

    summary_long = pd.DataFrame(coverage_summary_rows)
    summary_long, summary_wide = build_coverage_summary_tables(summary_long)

    summary_csv = os.path.join(run_outdir, "Source_Definer_Coverage_Summary.csv")
    summary_wide_csv = os.path.join(run_outdir, "Source_Definer_Coverage_Summary_Wide.csv")
    summary_xlsx = os.path.join(run_outdir, "Source_Definer_Coverage_Summary.xlsx")

    summary_long.to_csv(summary_csv, index=False)
    summary_wide.to_csv(summary_wide_csv, index=False)

    save_multi_sheet_xlsx(
        {
            "Long": summary_long,
            "Wide": summary_wide,
            "Detailed": detailed_df
        },
        summary_xlsx
    )

    return summary_long

# =============================================================================
# MAIN
# =============================================================================
def main():
    print("🚀 Starting within-clade source SNP analysis...")

    for p in [MASTER_REPORT, CLADE_META, SOURCE_FILE]:
        if not os.path.exists(p):
            print(f"❌ Missing required file: {p}")
            return

    df = pd.read_csv(MASTER_REPORT, low_memory=False)
    clade_meta = pd.read_csv(CLADE_META, sep="\t", low_memory=False)
    source_meta = parse_source_colorstrip(SOURCE_FILE)

    ensure_required_columns(df, ["Sample", "POS"])
    ensure_required_columns(clade_meta, ["Sample", "Lineage"])
    ensure_required_columns(source_meta, ["Sample", "Source"])

    df["Sample"] = df["Sample"].map(normalize_sample_id)
    clade_meta["Sample"] = clade_meta["Sample"].map(normalize_sample_id)
    source_meta["Sample"] = source_meta["Sample"].map(normalize_sample_id)

    clade_meta = clade_meta.dropna(subset=["Sample", "Lineage"]).drop_duplicates(subset=["Sample"]).copy()
    clade_meta = clade_meta[clade_meta["Lineage"].isin(FOCAL_CLADES)].copy()

    meta = clade_meta.merge(source_meta, on="Sample", how="inner")
    meta = meta.drop_duplicates(subset=["Sample"]).copy()

    print("\n📊 Metadata counts by clade and source:")
    print(meta.groupby(["Lineage", "Source"]).size())

    # ensure core columns exist
    for col in ["REF", "ALT", "GENE", "EFFECT", "AA_CHANGE"]:
        if col not in df.columns:
            df[col] = pd.NA

    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df = df.dropna(subset=["Sample", "POS"]).copy()
    df["POS"] = df["POS"].astype(int)

    # recombination flagging
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
        print("⚠️ Recombination GFF not found. Proceeding without recombination masking.")

    scenarios = [
        ("all_snps", df.copy()),
        ("vertical_only", df[df["is_recombination"] == False].copy())
    ]

    all_run_summaries = []

    for scenario_label, current_data in scenarios:
        run_outdir = os.path.join(OUTPUT_ROOT, scenario_label)
        print(f"\n🔵 Running scenario: {scenario_label}")
        run_summary = run_source_pipeline(meta, current_data, run_outdir)

        if run_summary is not None and not run_summary.empty:
            run_summary = run_summary.copy()
            run_summary["Scenario"] = scenario_label
            all_run_summaries.append(run_summary)

    # overall combined summaries
    if all_run_summaries:
        overall_long = pd.concat(all_run_summaries, ignore_index=True)
        overall_long["Coverage_Bin"] = pd.Categorical(
            overall_long["Coverage_Bin"],
            categories=COVERAGE_BIN_ORDER,
            ordered=True
        )
        overall_long = overall_long.sort_values(
            ["Scenario", "Clade", "Source_Target", "Coverage_Bin"]
        ).reset_index(drop=True)

        overall_csv = os.path.join(OUTPUT_ROOT, "Overall_Source_Definer_Coverage_Summary.csv")
        overall_wide_csv = os.path.join(OUTPUT_ROOT, "Overall_Source_Definer_Coverage_Summary_Wide.csv")
        overall_xlsx = os.path.join(OUTPUT_ROOT, "Overall_Source_Definer_Coverage_Summary.xlsx")

        overall_wide = overall_long.pivot_table(
            index=[
                "Scenario", "Clade", "Source_Target", "Samples_in_Clade",
                "Target_Samples", "Total_Significant_SNPs", "Total_Strict_Definers"
            ],
            columns="Coverage_Bin",
            values="Definer_SNP_Count",
            aggfunc="sum",
            fill_value=0
        ).reset_index()

        overall_wide.columns.name = None
        for c in COVERAGE_BIN_ORDER:
            if c not in overall_wide.columns:
                overall_wide[c] = 0

        overall_csv_by_clade_source = os.path.join(OUTPUT_ROOT, "Overall_Source_Definer_ByCladeSource.csv")
        overall_by_clade_source = overall_long.groupby(
            ["Clade", "Source_Target", "Coverage_Bin"],
            as_index=False
        ).agg(
            Definer_SNP_Count=("Definer_SNP_Count", "sum"),
            Runs_Contributing=("Scenario", "count")
        )

        overall_long.to_csv(overall_csv, index=False)
        overall_wide.to_csv(overall_wide_csv, index=False)
        overall_by_clade_source.to_csv(overall_csv_by_clade_source, index=False)

        save_multi_sheet_xlsx(
            {
                "By_Run_Long": overall_long,
                "By_Run_Wide": overall_wide,
                "By_CladeSource": overall_by_clade_source
            },
            overall_xlsx
        )

    print(f"\n🎉 ANALYSIS COMPLETE. Results located in: {OUTPUT_ROOT}")

if __name__ == "__main__":
    main()
