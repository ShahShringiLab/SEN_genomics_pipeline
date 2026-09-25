#!/usr/bin/env python3
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

# =============================================================================
# CONFIG
# =============================================================================
BASE = Path("/home/samuelajulo/SENBio/Final/iqtree_final/Snps_clade_integrated")
OUT_SUB = "master_defining_package"
TARGET_PREFIX = "Definer_Strict_"

# =============================================================================
# HELPERS
# =============================================================================
def norm(x) -> str:
    return "" if pd.isna(x) else str(x).strip()

def key(x) -> str:
    return re.sub(r"\s+", " ", norm(x).lower())

def is_blank(x) -> bool:
    return key(x) in {"", "nan", "none", ".", "-", "na", "n/a"}

def is_intergenic(x) -> bool:
    return key(x) == "intergenic"

def ensure_cols(df: pd.DataFrame, cols: List[str]) -> bool:
    return all(c in df.columns for c in cols)

def variant_id(r: pd.Series) -> str:
    return f"{norm(r.get('POS'))}|{norm(r.get('REF'))}|{norm(r.get('ALT'))}"

def wrap(seq: str, width: int = 80) -> str:
    seq = str(seq).strip()
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def classify_gene_annotation(x) -> str:
    s = norm(x)
    s_norm = key(x)
    if is_blank(s):
        return "Blank_or_missing"
    if s_norm == "intergenic":
        return "Intergenic"
    if s_norm == "hypothetical protein":
        return "Hypothetical_protein"
    if re.match(r"^SEN_RS\d+$", s):
        return "Locus_tag_like"
    return "Named_gene"

def infer_scheme_scenario_mode(csv_path: Path) -> Tuple[str, str, str]:
    """
    Expected path:
      .../Snps_clade_integrated/<scheme>/<scenario>/<mode>/Definer_Strict_<clade>.csv
    """
    parts = csv_path.parts
    idx = parts.index("Snps_clade_integrated")
    scheme = parts[idx + 1]
    scenario = parts[idx + 2]
    mode = parts[idx + 3]
    return scheme, scenario, mode

def infer_clade_from_stem(stem: str) -> str:
    if stem.startswith(TARGET_PREFIX):
        return stem[len(TARGET_PREFIX):].replace("_", " ")
    return stem.replace("_", " ")

# =============================================================================
# LOAD PROTEIN SOURCES
# =============================================================================
def load_reference_matches(path: Path) -> Dict[str, dict]:
    """
    Expected columns:
      query_gene, gene, locus_tag, old_locus_tag, protein_id, product, translation
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}

    df = pd.read_csv(path, sep="\t", low_memory=False)
    if df.empty or "query_gene" not in df.columns:
        return {}

    for c in ["gene", "locus_tag", "old_locus_tag", "protein_id", "product", "translation"]:
        if c not in df.columns:
            df[c] = pd.NA

    df["query_gene_norm"] = df["query_gene"].map(key)
    df["protein_sequence"] = df["translation"].fillna("").astype(str).str.strip()
    df["protein_status"] = df["protein_sequence"].apply(
        lambda x: "matched_reference" if x else "reference_missing_translation"
    )
    df["protein_length_aa"] = df["protein_sequence"].str.len()

    df = (
        df.sort_values(["query_gene_norm", "protein_length_aa"], ascending=[True, False])
          .drop_duplicates(subset=["query_gene_norm"], keep="first")
          .reset_index(drop=True)
    )

    return {row["query_gene_norm"]: row.to_dict() for _, row in df.iterrows()}

def build_rescue_status(row: pd.Series) -> str:
    flags = []
    if bool(row.get("partial_hit", False)):
        flags.append("partial_hit")
    if bool(row.get("internal_stop", False)):
        flags.append("internal_stop")
    if bool(row.get("frame_disruption", False)):
        flags.append("frame_disruption")
    if not flags:
        return "rescued_representative"
    return "rescued_representative;" + ";".join(flags)

def load_rescue_representatives(path: Path) -> Dict[str, dict]:
    """
    Expected columns from rescue representative output:
      query_gene, sample, sseqid, pident, qcov, bitscore,
      partial_hit, internal_stop, frame_disruption,
      longest_orf_aa, longest_orf_aa_len
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}

    df = pd.read_csv(path, sep="\t", low_memory=False)
    if df.empty or "query_gene" not in df.columns:
        return {}

    needed = [
        "sample", "sseqid", "pident", "qcov", "bitscore",
        "partial_hit", "internal_stop", "frame_disruption",
        "orf_found", "longest_orf_aa", "longest_orf_aa_len", "longest_orf_nt_len"
    ]
    for c in needed:
        if c not in df.columns:
            df[c] = pd.NA

    df["query_gene_norm"] = df["query_gene"].map(key)
    df["protein_sequence"] = df["longest_orf_aa"].fillna("").astype(str).str.strip()
    df["protein_length_aa"] = pd.to_numeric(df["longest_orf_aa_len"], errors="coerce").fillna(
        df["protein_sequence"].str.len()
    )
    df["protein_status"] = df.apply(build_rescue_status, axis=1)

    df = (
        df.sort_values(
            ["query_gene_norm", "protein_length_aa", "bitscore", "pident", "qcov"],
            ascending=[True, False, False, False, False]
        )
        .drop_duplicates(subset=["query_gene_norm"], keep="first")
        .reset_index(drop=True)
    )

    return {row["query_gene_norm"]: row.to_dict() for _, row in df.iterrows()}

# =============================================================================
# PROTEIN PICKER
# =============================================================================
def pick_protein(gene_value, ref_map: Dict[str, dict], rescue_map: Dict[str, dict]) -> dict:
    gene_raw = norm(gene_value)
    gene_norm = key(gene_value)

    out = {
        "Representative_Gene_Name": pd.NA,
        "Representative_Locus_Tag": pd.NA,
        "Representative_Old_Locus_Tag": pd.NA,
        "Representative_Protein_ID": pd.NA,
        "Representative_Product": pd.NA,
        "Protein_Source": "N/A",
        "Protein_Status": "no_gene_or_no_recoverable_protein",
        "Protein_Sequence": "N/A",
        "Protein_Length_AA": pd.NA,
        "Rescue_Sample": pd.NA,
        "Rescue_Contig": pd.NA,
        "Rescue_Pident": pd.NA,
        "Rescue_Qcov": pd.NA,
        "Rescue_Bitscore": pd.NA,
        "Rescue_Partial_Hit": pd.NA,
        "Rescue_Internal_Stop": pd.NA,
        "Rescue_Frame_Disruption": pd.NA,
    }

    if is_blank(gene_raw):
        out["Protein_Status"] = "blank_gene_annotation"
        return out

    if is_intergenic(gene_raw):
        out["Protein_Status"] = "intergenic_no_protein_expected"
        return out

    if gene_norm in ref_map:
        r = ref_map[gene_norm]
        seq = norm(r.get("protein_sequence"))
        if seq and seq != "N/A":
            out.update({
                "Representative_Gene_Name": r.get("gene", pd.NA),
                "Representative_Locus_Tag": r.get("locus_tag", pd.NA),
                "Representative_Old_Locus_Tag": r.get("old_locus_tag", pd.NA),
                "Representative_Protein_ID": r.get("protein_id", pd.NA),
                "Representative_Product": r.get("product", pd.NA),
                "Protein_Source": "reference",
                "Protein_Status": r.get("protein_status", "matched_reference"),
                "Protein_Sequence": seq,
                "Protein_Length_AA": len(seq),
            })
            return out

    if gene_norm in rescue_map:
        r = rescue_map[gene_norm]
        seq = norm(r.get("protein_sequence"))
        if seq and seq != "N/A":
            out.update({
                "Representative_Gene_Name": gene_raw,
                "Protein_Source": "rescued",
                "Protein_Status": r.get("protein_status", "rescued_representative"),
                "Protein_Sequence": seq,
                "Protein_Length_AA": len(seq),
                "Rescue_Sample": r.get("sample", pd.NA),
                "Rescue_Contig": r.get("sseqid", pd.NA),
                "Rescue_Pident": r.get("pident", pd.NA),
                "Rescue_Qcov": r.get("qcov", pd.NA),
                "Rescue_Bitscore": r.get("bitscore", pd.NA),
                "Rescue_Partial_Hit": r.get("partial_hit", pd.NA),
                "Rescue_Internal_Stop": r.get("internal_stop", pd.NA),
                "Rescue_Frame_Disruption": r.get("frame_disruption", pd.NA),
            })
            return out

    out["Protein_Status"] = "not_rescued_possible_premature_truncation_or_unresolved_cds"
    return out

# =============================================================================
# FASTA WRITER
# =============================================================================
def write_faa(df: pd.DataFrame, out_faa: Path) -> int:
    n = 0
    with open(out_faa, "w") as f:
        for _, r in df.iterrows():
            seq = norm(r.get("Protein_Sequence"))
            if not seq or seq == "N/A":
                continue

            header_parts = [
                f"gene_query={norm(r.get('GENE'))}",
                f"gene_key={norm(r.get('Gene_Key'))}",
                f"clade={norm(r.get('Clade'))}",
                f"scheme={norm(r.get('Scheme'))}",
                f"scenario={norm(r.get('Scenario'))}",
                f"mode={norm(r.get('Mode'))}",
                f"protein_source={norm(r.get('Protein_Source'))}",
                f"protein_status={norm(r.get('Protein_Status'))}",
                f"defining_snp_count={norm(r.get('Defining_SNP_Count_In_Gene'))}",
            ]

            for k in [
                "Representative_Locus_Tag",
                "Representative_Gene_Name",
                "Representative_Protein_ID",
                "Representative_Product",
            ]:
                val = norm(r.get(k))
                if val and val != "N/A":
                    header_parts.append(f"{k.lower()}={val}")

            f.write(">" + " | ".join(header_parts) + "\n")
            f.write(wrap(seq) + "\n")
            n += 1
    return n

# =============================================================================
# PROCESS ONE DEFINER FILE
# =============================================================================
def process_one_definer_file(fp: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(fp, low_memory=False)

    required = ["POS", "REF", "ALT", "GENE"]
    if not ensure_cols(df, required):
        raise ValueError(f"Raw definer file missing required columns {required}: {fp}")

    for c in ["EFFECT", "AA_CHANGE", "fdr", "percent_in", "percent_out", "Coverage_Bin"]:
        if c not in df.columns:
            df[c] = pd.NA

    scheme, scenario, mode = infer_scheme_scenario_mode(fp)
    clade = infer_clade_from_stem(fp.stem)

    protein_dir = fp.parent / "protein_sequence"
    rescue_dir = protein_dir / "query_rescue"

    ref_tsv = protein_dir / f"{fp.stem}_reference_matches.tsv"
    rescue_tsv = rescue_dir / f"{fp.stem}_query_rescue_representative.tsv"

    ref_map = load_reference_matches(ref_tsv)
    rescue_map = load_rescue_representatives(rescue_tsv)

    # SNP-level table
    snp_df = df.copy()
    snp_df["Scheme"] = scheme
    snp_df["Scenario"] = scenario
    snp_df["Mode"] = mode
    snp_df["Clade"] = clade
    snp_df["Source_File"] = str(fp)
    snp_df["Variant_ID"] = snp_df.apply(variant_id, axis=1)
    snp_df["Gene_Key"] = snp_df["GENE"].apply(lambda x: norm(x) if not is_blank(x) else "N/A")
    snp_df["Gene_Annotation_Class"] = snp_df["GENE"].apply(classify_gene_annotation)

    bundles = []
    for _, row in snp_df.iterrows():
        bundles.append(pick_protein(row.get("GENE"), ref_map, rescue_map))

    snp_df = pd.concat([snp_df.reset_index(drop=True), pd.DataFrame(bundles)], axis=1)

    snp_df["percent_in_num"] = pd.to_numeric(snp_df["percent_in"], errors="coerce")
    snp_df = snp_df.sort_values(
        ["Clade", "Gene_Key", "percent_in_num", "POS"],
        ascending=[True, True, False, True]
    ).drop(columns=["percent_in_num"]).reset_index(drop=True)

    # Gene-level table
    gene_rows = []
    for gene_key, g in snp_df.groupby("Gene_Key", dropna=False):
        first = g.iloc[0]
        gene_rows.append({
            "Scheme": first["Scheme"],
            "Scenario": first["Scenario"],
            "Mode": first["Mode"],
            "Clade": first["Clade"],
            "GENE": first["GENE"],
            "Gene_Key": gene_key,
            "Gene_Annotation_Class": first["Gene_Annotation_Class"],
            "Representative_Gene_Name": first["Representative_Gene_Name"],
            "Representative_Locus_Tag": first["Representative_Locus_Tag"],
            "Representative_Old_Locus_Tag": first["Representative_Old_Locus_Tag"],
            "Representative_Protein_ID": first["Representative_Protein_ID"],
            "Representative_Product": first["Representative_Product"],
            "Protein_Source": first["Protein_Source"],
            "Protein_Status": first["Protein_Status"],
            "Protein_Sequence": first["Protein_Sequence"],
            "Protein_Length_AA": first["Protein_Length_AA"],
            "Rescue_Sample": first["Rescue_Sample"],
            "Rescue_Contig": first["Rescue_Contig"],
            "Rescue_Pident": first["Rescue_Pident"],
            "Rescue_Qcov": first["Rescue_Qcov"],
            "Rescue_Bitscore": first["Rescue_Bitscore"],
            "Rescue_Partial_Hit": first["Rescue_Partial_Hit"],
            "Rescue_Internal_Stop": first["Rescue_Internal_Stop"],
            "Rescue_Frame_Disruption": first["Rescue_Frame_Disruption"],
            "Defining_SNP_Count_In_Gene": int(len(g)),
            "Has_Defining_SNP": 1,
            "Variant_IDs": ";".join(g["Variant_ID"].astype(str).tolist()),
            "Variant_Effects": ";".join(sorted([x for x in g["EFFECT"].fillna("").astype(str).unique() if x])),
            "Coverage_Bins": ";".join(sorted([x for x in g["Coverage_Bin"].fillna("").astype(str).unique() if x])),
            "Min_percent_in": pd.to_numeric(g["percent_in"], errors="coerce").min(),
            "Max_percent_in": pd.to_numeric(g["percent_in"], errors="coerce").max(),
            "Min_percent_out": pd.to_numeric(g["percent_out"], errors="coerce").min(),
            "Max_percent_out": pd.to_numeric(g["percent_out"], errors="coerce").max(),
        })

    gene_df = pd.DataFrame(gene_rows).sort_values(
        ["Clade", "Defining_SNP_Count_In_Gene", "Gene_Key"],
        ascending=[True, False, True]
    ).reset_index(drop=True)

    return snp_df, gene_df

# =============================================================================
# SAVE PER-ITERATION
# =============================================================================
def save_iteration_outputs(fp: Path, snp_df: pd.DataFrame, gene_df: pd.DataFrame):
    outdir = fp.parent / OUT_SUB
    outdir.mkdir(exist_ok=True)

    stem = fp.stem  # no .csv in middle anymore

    snp_tsv = outdir / f"{stem}_Master_Defining_SNPs.tsv"
    snp_csv = outdir / f"{stem}_Master_Defining_SNPs.csv"
    gene_tsv = outdir / f"{stem}_Master_Defining_Genes.tsv"
    gene_csv = outdir / f"{stem}_Master_Defining_Genes.csv"
    faa = outdir / f"{stem}_Master_Defining_Proteins.faa"

    snp_df.to_csv(snp_tsv, sep="\t", index=False)
    snp_df.to_csv(snp_csv, index=False)

    gene_df.to_csv(gene_tsv, sep="\t", index=False)
    gene_df.to_csv(gene_csv, index=False)

    n = write_faa(gene_df, faa)

    print(f"[OK] Iteration aggregated: {fp}")
    print(f"     SNP rows: {len(snp_df)}")
    print(f"     Gene rows: {len(gene_df)}")
    print(f"     FASTA proteins written: {n}")
    print(f"     -> {outdir}")

# =============================================================================
# SAVE FOLDER MASTER
# =============================================================================
def save_folder_master(folder: Path, all_snp: List[pd.DataFrame], all_gene: List[pd.DataFrame]):
    if not all_snp and not all_gene:
        return

    outdir = folder / OUT_SUB
    outdir.mkdir(exist_ok=True)

    folder_snp = pd.concat(all_snp, ignore_index=True) if all_snp else pd.DataFrame()
    folder_gene = pd.concat(all_gene, ignore_index=True) if all_gene else pd.DataFrame()

    snp_tsv = outdir / "Folder_Master_Defining_SNPs.tsv"
    snp_csv = outdir / "Folder_Master_Defining_SNPs.csv"
    gene_tsv = outdir / "Folder_Master_Defining_Genes.tsv"
    gene_csv = outdir / "Folder_Master_Defining_Genes.csv"
    faa = outdir / "Folder_Master_Defining_Proteins.faa"

    folder_snp.to_csv(snp_tsv, sep="\t", index=False)
    folder_snp.to_csv(snp_csv, index=False)

    folder_gene.to_csv(gene_tsv, sep="\t", index=False)
    folder_gene.to_csv(gene_csv, index=False)

    n = write_faa(folder_gene, faa)

    print(f"[OK] Folder master written: {folder}")
    print(f"     Folder SNP rows: {len(folder_snp)}")
    print(f"     Folder gene rows: {len(folder_gene)}")
    print(f"     Folder FASTA proteins written: {n}")

# =============================================================================
# SAVE OVERALL MASTER
# =============================================================================
def save_overall_master(base: Path, all_snp: List[pd.DataFrame], all_gene: List[pd.DataFrame]):
    outdir = base / OUT_SUB
    outdir.mkdir(exist_ok=True)

    overall_snp = pd.concat(all_snp, ignore_index=True) if all_snp else pd.DataFrame()
    overall_gene = pd.concat(all_gene, ignore_index=True) if all_gene else pd.DataFrame()

    snp_tsv = outdir / "Overall_Master_Defining_SNPs.tsv"
    snp_csv = outdir / "Overall_Master_Defining_SNPs.csv"
    gene_tsv = outdir / "Overall_Master_Defining_Genes.tsv"
    gene_csv = outdir / "Overall_Master_Defining_Genes.csv"
    faa = outdir / "Overall_Master_Defining_Proteins.faa"

    overall_snp.to_csv(snp_tsv, sep="\t", index=False)
    overall_snp.to_csv(snp_csv, index=False)

    overall_gene.to_csv(gene_tsv, sep="\t", index=False)
    overall_gene.to_csv(gene_csv, index=False)

    n = write_faa(overall_gene, faa)

    print("\n[OK] Overall master written")
    print(f"     Overall SNP rows: {len(overall_snp)}")
    print(f"     Overall gene rows: {len(overall_gene)}")
    print(f"     Overall FASTA proteins written: {n}")
    print(f"     -> {outdir}")

# =============================================================================
# FIND RAW DEFINER FILES ONLY
# =============================================================================
def find_raw_definer_files(base: Path) -> List[Path]:
    hits = []
    for root, dirs, files in os.walk(base):
        # do not descend into our own output folders
        dirs[:] = [d for d in dirs if d != OUT_SUB]

        for fn in files:
            if not fn.endswith(".csv"):
                continue
            if not fn.startswith(TARGET_PREFIX):
                continue
            if "_Master_" in fn:
                continue
            fp = Path(root) / fn

            # verify raw columns before accepting
            try:
                test_df = pd.read_csv(fp, nrows=3)
            except Exception:
                continue

            if ensure_cols(test_df, ["POS", "REF", "ALT", "GENE"]):
                hits.append(fp)

    return sorted(hits)

# =============================================================================
# MAIN
# =============================================================================
def main():
    if not BASE.exists():
        raise FileNotFoundError(f"Base folder not found: {BASE}")

    target_files = find_raw_definer_files(BASE)
    print(f"[INFO] Raw Definer_Strict files found: {len(target_files)}")

    if not target_files:
        print("[WARN] No valid raw Definer_Strict CSV files found.")
        return

    overall_snp = []
    overall_gene = []

    # group by parent folder
    files_by_folder: Dict[Path, List[Path]] = {}
    for fp in target_files:
        files_by_folder.setdefault(fp.parent, []).append(fp)

    for folder, fps in sorted(files_by_folder.items(), key=lambda x: str(x[0])):
        folder_snp = []
        folder_gene = []

        print(f"\n[INFO] Processing folder: {folder}")
        for fp in sorted(fps):
            try:
                snp_df, gene_df = process_one_definer_file(fp)
                save_iteration_outputs(fp, snp_df, gene_df)

                folder_snp.append(snp_df)
                folder_gene.append(gene_df)

                overall_snp.append(snp_df)
                overall_gene.append(gene_df)

            except Exception as e:
                print(f"[FAIL] {fp}: {e}")

        save_folder_master(folder, folder_snp, folder_gene)

    save_overall_master(BASE, overall_snp, overall_gene)

    print("\n✅ DONE")
    print("[INFO] Outputs created at three levels:")
    print("[INFO]   1. Per iteration (per Definer_Strict file)")
    print("[INFO]   2. Per folder (all clades combined in that folder)")
    print("[INFO]   3. Overall (all folders combined)")

if __name__ == "__main__":
    main()
