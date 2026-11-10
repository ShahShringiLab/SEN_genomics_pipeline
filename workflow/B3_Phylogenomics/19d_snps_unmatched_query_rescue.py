#!/usr/bin/env python3
import os
import re
import shutil
import tempfile
import subprocess
from pathlib import Path

import pandas as pd

# =============================================================================
# RESCUE UNMATCHED GENES ONLY USING CLADE-SPECIFIC ASSEMBLIES
#
# ADDED:
#   *_query_rescue_representative.faa
#   *_query_rescue_representative.tsv
#
# One best rescued protein per unmatched gene.
# =============================================================================

REPO_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = Path(os.environ.get("SEN_ROOT", REPO_ROOT))
REFERENCE_GBK = Path(os.environ.get("SEN_REFERENCE_GBK", BASE_DIR / "reference" / "reference.gbk"))
ASSEMBLY_DIR = Path(os.environ.get("SEN_FINAL_CONTIGS_DIR", BASE_DIR / "Final_Contigs_Only"))
SEARCH_ROOT = Path(os.environ.get("SEN_SNP_INTEGRATED_ROOT", BASE_DIR / "iqtree_final" / "Snps_clade_integrated"))
METADATA_FILE = Path(os.environ.get("SEN_CLADE_METADATA", BASE_DIR / "metadata" / "final_clade_metadata.tsv"))

THREADS_PER_BLAST = 1

MIN_IDENTITY = 90.0
MIN_QCOV = 80.0

UNMATCHED_SUFFIX = "_unmatched_genes.txt"

# =============================================================================
# HELPERS
# =============================================================================
def need_cmd(cmd: str):
    if shutil.which(cmd) is None:
        raise SystemExit(f"[ERROR] Required command not found in PATH: {cmd}")

def normalize_sample_name(x) -> str:
    s = str(x).strip()
    m = re.search(r"(SRR\d+)", s)
    if m:
        return m.group(1)
    return re.sub(r"\.(fasta|fa|fna)$", "", os.path.basename(s), flags=re.IGNORECASE)

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def translate_dna(seq: str) -> str:
    seq = seq.upper()
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)

def longest_orf_from_frame(seq: str):
    seq = seq.upper()
    best = {
        "start": None,
        "end": None,
        "nt_seq": "",
        "aa_seq": "",
        "length_nt": 0,
        "frame": None,
        "internal_stops": None,
    }

    stop_codons = {"TAA", "TAG", "TGA"}

    for frame in [0, 1, 2]:
        i = frame
        while i <= len(seq) - 3:
            codon = seq[i:i+3]
            if codon == "ATG":
                j = i
                found_stop = False
                while j <= len(seq) - 3:
                    c = seq[j:j+3]
                    if c in stop_codons:
                        nt = seq[i:j+3]
                        aa = translate_dna(nt)
                        if len(nt) > best["length_nt"]:
                            best = {
                                "start": i + 1,
                                "end": j + 3,
                                "nt_seq": nt,
                                "aa_seq": aa,
                                "length_nt": len(nt),
                                "frame": frame + 1,
                                "internal_stops": aa[:-1].count("*"),
                            }
                        found_stop = True
                        break
                    j += 3
                i = j + 3 if found_stop else i + 3
            else:
                i += 3

    return best

def read_fasta_dict(path: Path):
    seqs = {}
    header = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            seqs[header] = "".join(chunks)
    return seqs

def parse_reference_cds_from_gbk(gbk_path: Path):
    with open(gbk_path, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    genome_chunks = []
    in_origin = False
    for line in lines:
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if line.startswith("//"):
                break
            seq_part = re.sub(r"[^acgtACGTnN]", "", line)
            genome_chunks.append(seq_part)
    genome_seq = "".join(genome_chunks).upper()

    records = []
    current = None
    in_cds = False

    def flush():
        nonlocal current
        if current is not None:
            if current.get("coords") and (
                current.get("locus_tag") or current.get("gene") or current.get("protein_id")
            ):
                try:
                    nt_seq = extract_cds_from_coords(genome_seq, current["coords"])
                except Exception:
                    nt_seq = None
                current["cds_nt"] = nt_seq
                records.append(current)
        current = None

    for i, raw in enumerate(lines):
        line = raw.rstrip("\n")
        if line.startswith("     CDS"):
            flush()
            in_cds = True
            coords = line[21:].strip()
            current = {
                "coords": coords,
                "gene": None,
                "locus_tag": None,
                "old_locus_tag": None,
                "protein_id": None,
                "product": None,
            }
            continue

        if in_cds and re.match(r"^     [A-Za-z_]", line) and not line.startswith("                     /"):
            flush()
            in_cds = False
            continue

        if in_cds and line.startswith("                     /"):
            text = line.strip()

            def collect_value(start_text, start_idx):
                if '="' not in start_text:
                    return None, start_idx
                value = start_text.split('="', 1)[1]
                if value.endswith('"'):
                    return value[:-1], start_idx
                parts = [value]
                j = start_idx + 1
                while j < len(lines):
                    nxt = lines[j].rstrip("\n")
                    if nxt.startswith("                     "):
                        piece = nxt.strip()
                        if piece.endswith('"'):
                            parts.append(piece[:-1])
                            return "".join(parts), j
                        parts.append(piece)
                        j += 1
                    else:
                        return "".join(parts), j - 1
                return "".join(parts), j - 1

            if text.startswith("/gene="):
                current["gene"], _ = collect_value(text, i)
            elif text.startswith("/locus_tag="):
                current["locus_tag"], _ = collect_value(text, i)
            elif text.startswith("/old_locus_tag="):
                current["old_locus_tag"], _ = collect_value(text, i)
            elif text.startswith("/protein_id="):
                current["protein_id"], _ = collect_value(text, i)
            elif text.startswith("/product="):
                current["product"], _ = collect_value(text, i)

    flush()

    df = pd.DataFrame(records)
    if df.empty:
        raise RuntimeError(f"No CDS records parsed from {gbk_path}")

    for col in ["gene", "locus_tag", "old_locus_tag", "protein_id"]:
        if col not in df.columns:
            df[col] = None
        df[f"{col}_norm"] = df[col].fillna("").astype(str).str.strip().str.lower()

    return df

def extract_cds_from_coords(genome_seq: str, coord_text: str) -> str:
    s = coord_text.replace("<", "").replace(">", "").strip()
    strand = 1

    if s.startswith("complement(") and s.endswith(")"):
        strand = -1
        s = s[len("complement("):-1].strip()

    if s.startswith("join(") and s.endswith(")"):
        s = s[len("join("):-1].strip()
        parts = [p.strip() for p in s.split(",")]
    else:
        parts = [s]

    seq_chunks = []
    for part in parts:
        m = re.match(r"(\d+)\.\.(\d+)", part)
        if not m:
            continue
        start = int(m.group(1))
        end = int(m.group(2))
        seq_chunks.append(genome_seq[start-1:end])

    seq = "".join(seq_chunks)
    if strand == -1:
        seq = revcomp(seq)
    return seq

def find_unmatched_files(search_root: Path):
    return sorted(search_root.rglob(f"*{UNMATCHED_SUFFIX}"))

def read_unmatched_genes(path: Path):
    genes = []
    with open(path) as f:
        for line in f:
            g = line.strip()
            if not g:
                continue
            if g.lower() in {"intergenic", "nan", "none", ".", "-"}:
                continue
            genes.append(g)
    return sorted(set(genes))

def match_ref_cds(genes, ref_df):
    matched = []
    unmatched = []
    for g in genes:
        q = g.strip().lower()
        hit = ref_df[ref_df["locus_tag_norm"] == q]
        if hit.empty:
            hit = ref_df[ref_df["old_locus_tag_norm"] == q]
        if hit.empty:
            hit = ref_df[ref_df["gene_norm"] == q]
        if hit.empty:
            hit = ref_df[ref_df["protein_id_norm"] == q]

        if hit.empty:
            unmatched.append(g)
            continue

        row = hit.iloc[0]
        if pd.isna(row.get("cds_nt")) or not str(row.get("cds_nt")).strip():
            unmatched.append(g)
            continue

        matched.append({
            "query_gene": g,
            "gene": row.get("gene"),
            "locus_tag": row.get("locus_tag"),
            "old_locus_tag": row.get("old_locus_tag"),
            "protein_id": row.get("protein_id"),
            "product": row.get("product"),
            "cds_nt": str(row.get("cds_nt")).upper(),
            "ref_cds_len": len(str(row.get("cds_nt")).upper()),
        })
    return matched, unmatched

def write_query_fasta(ref_records, out_fasta):
    with open(out_fasta, "w") as f:
        for r in ref_records:
            header = r["query_gene"].replace(" ", "_")
            f.write(f">{header}\n")
            seq = r["cds_nt"]
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

def build_blast_db(genome_path: Path, db_prefix: Path):
    cmd = [
        "makeblastdb",
        "-in", str(genome_path),
        "-dbtype", "nucl",
        "-out", str(db_prefix),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def run_blast(query_fasta: Path, db_prefix: Path, out_tsv: Path):
    cmd = [
        "blastn",
        "-task", "blastn",
        "-query", str(query_fasta),
        "-db", str(db_prefix),
        "-out", str(out_tsv),
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen",
        "-num_threads", str(THREADS_PER_BLAST),
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def extract_subject_region(genome_fasta: Path, contig_id: str, sstart: int, send: int):
    seqs = read_fasta_dict(genome_fasta)
    if contig_id not in seqs:
        return None, None
    contig_seq = seqs[contig_id]
    start = min(sstart, send)
    end = max(sstart, send)
    region = contig_seq[start-1:end]
    strand = "+" if sstart <= send else "-"
    if strand == "-":
        region = revcomp(region)
    return region.upper(), strand

def analyze_hit_orf(hit_nt: str, ref_cds_nt: str):
    result = {
        "partial_hit": False,
        "internal_stop": False,
        "frame_disruption": False,
        "orf_found": False,
        "longest_orf_nt_len": 0,
        "longest_orf_aa_len": 0,
        "longest_orf_frame": None,
        "longest_orf_nt": "",
        "longest_orf_aa": "",
    }

    if not hit_nt:
        result["partial_hit"] = True
        result["frame_disruption"] = True
        return result

    if len(hit_nt) < len(ref_cds_nt):
        result["partial_hit"] = True

    orf = longest_orf_from_frame(hit_nt)
    if orf["length_nt"] == 0:
        result["frame_disruption"] = True
        return result

    result["orf_found"] = True
    result["longest_orf_nt_len"] = orf["length_nt"]
    result["longest_orf_aa_len"] = len(orf["aa_seq"])
    result["longest_orf_frame"] = orf["frame"]
    result["longest_orf_nt"] = orf["nt_seq"]
    result["longest_orf_aa"] = orf["aa_seq"]

    if orf["internal_stops"] and orf["internal_stops"] > 0:
        result["internal_stop"] = True

    if orf["length_nt"] % 3 != 0:
        result["frame_disruption"] = True

    if orf["length_nt"] < 0.8 * len(ref_cds_nt):
        result["frame_disruption"] = True

    return result

def write_faa(records, out_faa):
    n = 0
    with open(out_faa, "w") as f:
        for r in records:
            seq = str(r.get("longest_orf_aa", "")).strip()
            if not seq:
                continue
            header = [
                f"query_gene={r['query_gene']}",
                f"sample={r['sample']}",
                f"contig={r['sseqid']}",
                f"pident={r['pident']}",
                f"qcov={r['qcov']}",
                f"partial_hit={r['partial_hit']}",
                f"internal_stop={r['internal_stop']}",
                f"frame_disruption={r['frame_disruption']}",
            ]
            f.write(">" + " | ".join(map(str, header)) + "\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")
            n += 1
    return n

# =============================================================================
# CLADE-SPECIFIC ASSEMBLY SELECTION
# =============================================================================
def load_metadata(meta_file: Path):
    meta = pd.read_csv(meta_file, sep="\t", low_memory=False)
    if "Sample" not in meta.columns or "Lineage" not in meta.columns:
        raise SystemExit("[ERROR] Clade_metadata.txt must contain Sample and Lineage columns")
    meta["Sample"] = meta["Sample"].map(normalize_sample_name)
    meta = meta.dropna(subset=["Sample", "Lineage"]).copy()
    meta = meta[meta["Lineage"] != "Reference"].copy()
    return meta

def map_assemblies_by_sample(assembly_dir: Path):
    assembly_map = {}
    for ext in ("*.fasta", "*.fa", "*.fna"):
        for fp in assembly_dir.glob(ext):
            sample = normalize_sample_name(fp)
            assembly_map[sample] = fp
    return assembly_map

def infer_target_clades(unmatched_file: Path):
    path_str = str(unmatched_file)
    is_combined = "Combined_Clade1" in path_str
    is_split = "Split_Clade1A1B" in path_str
    fname = unmatched_file.name.replace("_unmatched_genes.txt", "")

    if "Clade_1A" in fname:
        return ["Clade 1A"]
    if "Clade_1B" in fname:
        return ["Clade 1B"]
    if "Clade_2" in fname:
        return ["Clade 2"]
    if "Clade_3" in fname:
        return ["Clade 3"]
    if "Unassigned" in fname:
        return ["Unassigned"]

    if is_combined and "Clade_1" in fname:
        return ["Clade 1A", "Clade 1B"]

    if is_split and "Clade_1" in fname:
        return ["Clade 1A", "Clade 1B"]

    return None

def select_clade_specific_assemblies(unmatched_file: Path, meta_df: pd.DataFrame, assembly_map: dict):
    target_clades = infer_target_clades(unmatched_file)
    if target_clades is None:
        return [], []

    samples = meta_df.loc[meta_df["Lineage"].isin(target_clades), "Sample"].dropna().unique().tolist()
    assemblies = [assembly_map[s] for s in samples if s in assembly_map]

    return target_clades, sorted(assemblies, key=lambda x: x.name)

# =============================================================================
# REPRESENTATIVE SELECTION
# =============================================================================
def select_representative_hits(hits_df: pd.DataFrame) -> pd.DataFrame:
    if hits_df.empty:
        return hits_df.copy()

    rep = hits_df.copy()

    # prefer biologically cleaner sequences first
    rep["_frame_ok"] = (~rep["frame_disruption"].fillna(True)).astype(int)
    rep["_stop_ok"] = (~rep["internal_stop"].fillna(True)).astype(int)
    rep["_partial_ok"] = (~rep["partial_hit"].fillna(True)).astype(int)
    rep["_orf_ok"] = rep["orf_found"].fillna(False).astype(int)

    rep = rep.sort_values(
        by=[
            "query_gene",
            "_frame_ok",
            "_stop_ok",
            "_partial_ok",
            "_orf_ok",
            "bitscore",
            "pident",
            "qcov",
            "longest_orf_aa_len",
            "longest_orf_nt_len",
            "sample",
        ],
        ascending=[True, False, False, False, False, False, False, False, False, False, True]
    )

    rep = rep.drop_duplicates(subset=["query_gene"], keep="first").copy()

    rep = rep.drop(columns=["_frame_ok", "_stop_ok", "_partial_ok", "_orf_ok"])
    return rep

# =============================================================================
# MAIN PROCESSOR
# =============================================================================
def process_one_unmatched_file(unmatched_file: Path, ref_df: pd.DataFrame, meta_df: pd.DataFrame, assembly_map: dict):
    genes = read_unmatched_genes(unmatched_file)
    if not genes:
        print(f"[WARN] No usable unmatched genes in {unmatched_file}")
        return

    target_clades, assemblies = select_clade_specific_assemblies(unmatched_file, meta_df, assembly_map)
    if not target_clades:
        print(f"[WARN] Could not infer target clade(s) from {unmatched_file}")
        return
    if not assemblies:
        print(f"[WARN] No clade-specific assemblies found for {unmatched_file} | target_clades={target_clades}")
        return

    ref_records, still_unmatched = match_ref_cds(genes, ref_df)

    out_dir = unmatched_file.parent
    stem = unmatched_file.name.replace("_unmatched_genes.txt", "")
    rescue_dir = out_dir / "query_rescue"
    rescue_dir.mkdir(exist_ok=True)

    ref_query_fasta = rescue_dir / f"{stem}_unmatched_reference_cds.fna"
    ref_match_table = rescue_dir / f"{stem}_unmatched_reference_cds_matches.tsv"
    still_unmatched_txt = rescue_dir / f"{stem}_still_unmatched_after_reference_cds.txt"
    rescue_hits_tsv = rescue_dir / f"{stem}_query_rescue_hits.tsv"
    rescue_faa = rescue_dir / f"{stem}_query_rescue.faa"
    rescue_meta_txt = rescue_dir / f"{stem}_query_rescue_search_scope.txt"

    # NEW representative outputs
    rescue_rep_tsv = rescue_dir / f"{stem}_query_rescue_representative.tsv"
    rescue_rep_faa = rescue_dir / f"{stem}_query_rescue_representative.faa"

    pd.DataFrame(ref_records).to_csv(ref_match_table, sep="\t", index=False)
    with open(still_unmatched_txt, "w") as f:
        for g in still_unmatched:
            f.write(g + "\n")

    with open(rescue_meta_txt, "w") as f:
        f.write(f"Target file: {unmatched_file}\n")
        f.write(f"Target clades: {', '.join(target_clades)}\n")
        f.write(f"Assemblies searched: {len(assemblies)}\n")
        for a in assemblies:
            f.write(f"{normalize_sample_name(a)}\t{a}\n")

    if not ref_records:
        print(f"[WARN] No reference CDS sequences recovered for {unmatched_file}")
        return

    write_query_fasta(ref_records, ref_query_fasta)
    ref_map = {r["query_gene"]: r for r in ref_records}

    all_rows = []

    with tempfile.TemporaryDirectory(prefix="blast_rescue_") as tmpdir:
        tmpdir = Path(tmpdir)

        for genome in assemblies:
            sample = normalize_sample_name(genome)
            db_prefix = tmpdir / f"{sample}_db"
            out_tsv = tmpdir / f"{sample}.tsv"

            try:
                build_blast_db(genome, db_prefix)
                run_blast(ref_query_fasta, db_prefix, out_tsv)

                if not out_tsv.exists() or out_tsv.stat().st_size == 0:
                    continue

                df = pd.read_csv(
                    out_tsv,
                    sep="\t",
                    header=None,
                    names=[
                        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                        "qlen", "slen"
                    ]
                )
                if df.empty:
                    continue

                df["sample"] = sample
                df["qcov"] = 100.0 * df["length"] / df["qlen"]
                df = df[(df["pident"] >= MIN_IDENTITY) & (df["qcov"] >= MIN_QCOV)].copy()
                if df.empty:
                    continue

                # top hit per sample per gene
                df = df.sort_values(
                    ["sample", "qseqid", "bitscore", "pident", "qcov", "length"],
                    ascending=[True, True, False, False, False, False]
                ).drop_duplicates(subset=["sample", "qseqid"], keep="first")

                for _, row in df.iterrows():
                    hit_nt, strand = extract_subject_region(
                        genome,
                        row["sseqid"],
                        int(row["sstart"]),
                        int(row["send"])
                    )
                    ref_nt = ref_map[row["qseqid"]]["cds_nt"]
                    orf_info = analyze_hit_orf(hit_nt, ref_nt)

                    out_row = {
                        "sample": row["sample"],
                        "query_gene": row["qseqid"],
                        "sseqid": row["sseqid"],
                        "pident": float(row["pident"]),
                        "qcov": float(row["qcov"]),
                        "bitscore": float(row["bitscore"]),
                        "evalue": float(row["evalue"]),
                        "qstart": int(row["qstart"]),
                        "qend": int(row["qend"]),
                        "sstart": int(row["sstart"]),
                        "send": int(row["send"]),
                        "strand": strand,
                        "hit_nt_len": len(hit_nt) if hit_nt else 0,
                        "ref_cds_len": len(ref_nt),
                        **orf_info
                    }
                    all_rows.append(out_row)

            except Exception as e:
                print(f"[WARN] Failed rescue search for {sample} in {unmatched_file.name}: {e}")

    hits_df = pd.DataFrame(all_rows)
    if hits_df.empty:
        print(f"[WARN] No rescued hits passing filters for {unmatched_file}")
        return

    hits_df = hits_df.sort_values(["query_gene", "sample", "bitscore"], ascending=[True, True, False]).reset_index(drop=True)
    hits_df.to_csv(rescue_hits_tsv, sep="\t", index=False)

    faa_records = hits_df.to_dict("records")
    n_faa = write_faa(faa_records, rescue_faa)

    # NEW representative outputs
    rep_df = select_representative_hits(hits_df)
    rep_df.to_csv(rescue_rep_tsv, sep="\t", index=False)
    n_rep_faa = write_faa(rep_df.to_dict("records"), rescue_rep_faa)

    print(f"[OK] {unmatched_file}")
    print(f"     target clades: {target_clades}")
    print(f"     clade-specific assemblies searched: {len(assemblies)}")
    print(f"     unmatched genes listed: {len(genes)}")
    print(f"     reference CDS recovered for search: {len(ref_records)}")
    print(f"     still unmatched after ref CDS lookup: {len(still_unmatched)}")
    print(f"     rescued top hits: {len(hits_df)}")
    print(f"     rescued faa entries: {n_faa}")
    print(f"     representative rescued genes: {len(rep_df)}")
    print(f"     representative faa entries: {n_rep_faa}")
    print(f"     -> {rescue_dir}")

def main():
    need_cmd("makeblastdb")
    need_cmd("blastn")

    if not REFERENCE_GBK.exists():
        raise FileNotFoundError(f"Reference GBK not found: {REFERENCE_GBK}")
    if not ASSEMBLY_DIR.exists():
        raise FileNotFoundError(f"Assembly directory not found: {ASSEMBLY_DIR}")
    if not SEARCH_ROOT.exists():
        raise FileNotFoundError(f"Search root not found: {SEARCH_ROOT}")
    if not METADATA_FILE.exists():
        raise FileNotFoundError(f"Metadata file not found: {METADATA_FILE}")

    ref_df = parse_reference_cds_from_gbk(REFERENCE_GBK)
    print(f"[INFO] Reference CDS parsed: {len(ref_df)}")

    meta_df = load_metadata(METADATA_FILE)
    assembly_map = map_assemblies_by_sample(ASSEMBLY_DIR)
    print(f"[INFO] Assemblies indexed by sample: {len(assembly_map)}")

    unmatched_files = find_unmatched_files(SEARCH_ROOT)
    print(f"[INFO] Unmatched gene files found: {len(unmatched_files)}")

    if not unmatched_files:
        print("[WARN] No *_unmatched_genes.txt files found.")
        return

    for fp in unmatched_files:
        process_one_unmatched_file(fp, ref_df, meta_df, assembly_map)

    print("\n✅ DONE")
    print("[INFO] Existing matched reference outputs were not changed.")
    print("[INFO] Rescue outputs are in each protein_sequence/query_rescue/ folder.")
    print("[INFO] Representative rescue files were also created:")
    print("[INFO]   *_query_rescue_representative.tsv")
    print("[INFO]   *_query_rescue_representative.faa")

if __name__ == "__main__":
    main()
