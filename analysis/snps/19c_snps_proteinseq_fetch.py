#!/usr/bin/env python3
import os
import re
from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SEN_ROOT = Path(os.environ.get("SEN_ROOT", REPO_ROOT))
GBK_FILE = os.environ.get("SEN_REFERENCE_GBK", str(SEN_ROOT / "reference" / "reference.gbk"))
SEARCH_ROOT = os.environ.get("SEN_SNP_INTEGRATED_ROOT", str(SEN_ROOT / "iqtree_final" / "Snps_clade_integrated"))
TARGET_PREFIXES = ("Definer_Strict_", "Significant_All_")

def parse_genbank_proteins(gbk_path):
    """
    Robust parser for NCBI/RefSeq GenBank CDS features with /translation qualifiers.
    Extracts:
      gene, locus_tag, old_locus_tag, protein_id, product, translation
    """
    records = []
    current = None
    in_cds = False

    def flush():
        nonlocal current
        if current and current.get("translation"):
            records.append(current)
        current = None

    with open(gbk_path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].rstrip("\n")

        # start of a CDS feature
        if line.startswith("     CDS"):
            flush()
            in_cds = True
            current = {
                "gene": None,
                "locus_tag": None,
                "old_locus_tag": None,
                "protein_id": None,
                "product": None,
                "translation": None,
            }
            i += 1
            continue

        # new top-level feature starts, so CDS block ends
        if in_cds and re.match(r"^     [A-Za-z_]", line) and not line.startswith("                     /"):
            flush()
            in_cds = False
            continue

        if in_cds and line.startswith("                     /"):
            text = line.strip()

            def collect_multiline_quoted(start_text, start_index):
                """
                Collect qualifier values that may span multiple lines.
                start_text begins with /key="...
                """
                if '="' not in start_text:
                    return None, start_index

                value = start_text.split('="', 1)[1]

                # one-line qualifier
                if value.endswith('"'):
                    return value[:-1], start_index

                chunks = [value]
                j = start_index + 1
                while j < len(lines):
                    nxt = lines[j].rstrip("\n")

                    # continuation lines in GenBank qualifiers are indented 21 spaces
                    if nxt.startswith("                     "):
                        piece = nxt.strip()
                        if piece.endswith('"'):
                            chunks.append(piece[:-1])
                            return "".join(chunks), j
                        else:
                            chunks.append(piece)
                            j += 1
                            continue
                    else:
                        # malformed or truncated qualifier
                        return "".join(chunks), j - 1

                return "".join(chunks), j - 1

            if text.startswith("/gene="):
                val, i = collect_multiline_quoted(text, i)
                current["gene"] = val

            elif text.startswith("/locus_tag="):
                val, i = collect_multiline_quoted(text, i)
                current["locus_tag"] = val

            elif text.startswith("/old_locus_tag="):
                val, i = collect_multiline_quoted(text, i)
                current["old_locus_tag"] = val

            elif text.startswith("/protein_id="):
                val, i = collect_multiline_quoted(text, i)
                current["protein_id"] = val

            elif text.startswith("/product="):
                val, i = collect_multiline_quoted(text, i)
                current["product"] = val

            elif text.startswith("/translation="):
                val, i = collect_multiline_quoted(text, i)
                if val is not None:
                    val = re.sub(r"\s+", "", val).upper()
                    current["translation"] = val if val else None

        i += 1

    flush()

    df = pd.DataFrame(records)
    if df.empty:
        raise RuntimeError(f"No CDS translations found in {gbk_path}")

    for col in ["gene", "locus_tag", "old_locus_tag", "protein_id"]:
        if col not in df.columns:
            df[col] = None
        df[f"{col}_norm"] = df[col].fillna("").astype(str).str.strip().str.lower()

    return df

def clean_gene_list(df):
    if "GENE" not in df.columns:
        return []

    genes = (
        df["GENE"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )

    bad = {"", "nan", "none", ".", "-", "intergenic", "hypothetical protein"}
    genes = [g for g in genes if g.lower() not in bad]
    return sorted(set(genes))

def match_reference_proteins(genes, ref_df):
    """
    Match priority:
      1. locus_tag
      2. old_locus_tag
      3. gene
      4. protein_id
    """
    matched = []
    unmatched = []
    seen = set()

    for g in genes:
        g_norm = g.strip().lower()

        hit = ref_df[ref_df["locus_tag_norm"] == g_norm]
        if hit.empty:
            hit = ref_df[ref_df["old_locus_tag_norm"] == g_norm]
        if hit.empty:
            hit = ref_df[ref_df["gene_norm"] == g_norm]
        if hit.empty:
            hit = ref_df[ref_df["protein_id_norm"] == g_norm]

        if hit.empty:
            unmatched.append(g)
            continue

        row = hit.iloc[0]
        uniq = row["locus_tag"] if pd.notna(row["locus_tag"]) else (row["gene"] if pd.notna(row["gene"]) else row["protein_id"])
        if uniq in seen:
            continue
        seen.add(uniq)

        matched.append({
            "query_gene": g,
            "gene": row.get("gene"),
            "locus_tag": row.get("locus_tag"),
            "old_locus_tag": row.get("old_locus_tag"),
            "protein_id": row.get("protein_id"),
            "product": row.get("product"),
            "translation": row.get("translation"),
        })

    return matched, unmatched

def write_faa(records, out_faa):
    n = 0
    with open(out_faa, "w") as f:
        for r in records:
            seq = str(r["translation"]).strip()
            if not seq or seq.lower() == "none":
                continue

            header_parts = []
            for key in ["locus_tag", "old_locus_tag", "gene", "protein_id", "product"]:
                val = r.get(key)
                if pd.notna(val) and str(val).strip():
                    header_parts.append(f"{key}={str(val).strip()}")
            header_parts.append(f"query={r['query_gene']}")

            f.write(">" + " | ".join(header_parts) + "\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
            n += 1
    return n

def process_csv(csv_path, ref_df):
    df = pd.read_csv(csv_path)
    genes = clean_gene_list(df)

    out_dir = os.path.join(os.path.dirname(csv_path), "protein_sequence")
    os.makedirs(out_dir, exist_ok=True)

    stem = os.path.splitext(os.path.basename(csv_path))[0]
    out_faa = os.path.join(out_dir, f"{stem}_reference.faa")
    out_unmatched = os.path.join(out_dir, f"{stem}_unmatched_genes.txt")
    out_match_table = os.path.join(out_dir, f"{stem}_reference_matches.tsv")

    matched, unmatched = match_reference_proteins(genes, ref_df)
    n_written = write_faa(matched, out_faa)

    pd.DataFrame(matched).to_csv(out_match_table, sep="\t", index=False)

    with open(out_unmatched, "w") as f:
        for g in unmatched:
            f.write(g + "\n")

    print(f"[OK] {csv_path}")
    print(f"     genes in csv: {len(genes)}")
    print(f"     matched ref proteins: {len(matched)}")
    print(f"     faa written: {n_written}")
    print(f"     unmatched: {len(unmatched)}")
    print(f"     -> {out_faa}")

def main():
    if not os.path.exists(GBK_FILE):
        raise FileNotFoundError(f"Reference GBK not found: {GBK_FILE}")

    ref_df = parse_genbank_proteins(GBK_FILE)
    print(f"[INFO] Reference CDS proteins parsed: {len(ref_df)}")

    target_csvs = []
    for root, _, files in os.walk(SEARCH_ROOT):
        for fn in files:
            if fn.endswith(".csv") and fn.startswith(TARGET_PREFIXES):
                target_csvs.append(os.path.join(root, fn))

    target_csvs = sorted(target_csvs)

    if not target_csvs:
        print("[WARN] No target CSV files found.")
        return

    print(f"[INFO] Target CSV files found: {len(target_csvs)}")

    for csv_path in target_csvs:
        try:
            process_csv(csv_path, ref_df)
        except Exception as e:
            print(f"[FAIL] {csv_path}: {e}")

if __name__ == "__main__":
    main()
