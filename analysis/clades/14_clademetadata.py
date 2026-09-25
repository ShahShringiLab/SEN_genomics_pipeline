import os
import re
import pandas as pd

# =============================================================================
# Updated behavior:
#  - Extract ONLY true leaf IDs from .treefile (SRR... + optional REF_NAME)
#  - Avoids bootstrap values contaminating leaf extraction
#  - Reads ITOL clade lists:
#       Clade1.txt
#       Clade1b.txt
#       Clade2.txt
#       Clade3.txt
#  - Creates:
#       Clade 1A = Clade1 - Clade1b
#  - Writes:
#       1) Clade_metadata.txt
#       2) Unassigned.txt
#       3) Clade_counts.txt
# =============================================================================

# 1) Paths
BASE_DIR = "/home/samuelajulo/SENBio/Final"
ITOL_DIR = os.path.join(BASE_DIR, "ITOL")
TREEFILE = os.path.join(BASE_DIR, "iqtree_final", "SSLAB_FINAL.treefile")

# Input clade files
CLADE1_FILE = os.path.join(ITOL_DIR, "Clade1.txt")
CLADE1B_FILE = os.path.join(ITOL_DIR, "Clade1b.txt")
CLADE2_FILE = os.path.join(ITOL_DIR, "Clade2.txt")
CLADE3_FILE = os.path.join(ITOL_DIR, "Clade3.txt")

# Outputs
OUT_METADATA = os.path.join(ITOL_DIR, "Clade_metadata.txt")
OUT_UNASSIGNED = os.path.join(ITOL_DIR, "Unassigned.txt")
OUT_COUNTS = os.path.join(ITOL_DIR, "Clade_counts.txt")

# Optional reference label
REF_NAME = "Reference"  # set to None if not needed


# -----------------------------------------------------------------------------
def load_srr_list(path: str) -> set:
    """Load SRR IDs from a text file, robustly."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        raise FileNotFoundError(f"Missing/empty: {path}")
    with open(path, "r") as f:
        txt = f.read()
    return set(re.findall(r"SRR\d+", txt))


def extract_leaves_from_tree(tree_path: str, ref_name=None) -> set:
    """Extract leaf IDs: SRR IDs + optional reference token."""
    if not os.path.exists(tree_path) or os.path.getsize(tree_path) == 0:
        raise FileNotFoundError(f"Missing/empty treefile: {tree_path}")

    with open(tree_path, "r") as f:
        tree = f.read()

    leaves = set(re.findall(r"SRR\d+", tree))

    if ref_name:
        if re.search(rf"(?<![A-Za-z0-9_.-]){re.escape(ref_name)}(?![A-Za-z0-9_.-])", tree):
            leaves.add(ref_name)

    return leaves


# -----------------------------------------------------------------------------
# 2) Extract all leaf IDs safely
all_leaves = extract_leaves_from_tree(TREEFILE, REF_NAME)

# 3) Load clades
clade1 = load_srr_list(CLADE1_FILE)
clade1b = load_srr_list(CLADE1B_FILE)
clade2 = load_srr_list(CLADE2_FILE)
clade3 = load_srr_list(CLADE3_FILE)

# 4) Define Clade 1A
clade1a = clade1 - clade1b

# 5) Sanity checks
missing_1b_from_1 = sorted(clade1b - clade1)
if missing_1b_from_1:
    print(f"[WARN] {len(missing_1b_from_1)} IDs in Clade1b are not present in Clade1.")
    print("       They will still be labeled as Clade 1B.")

overlaps = {
    "Clade1A_vs_Clade1B": clade1a & clade1b,
    "Clade1A_vs_Clade2": clade1a & clade2,
    "Clade1A_vs_Clade3": clade1a & clade3,
    "Clade1B_vs_Clade2": clade1b & clade2,
    "Clade1B_vs_Clade3": clade1b & clade3,
    "Clade2_vs_Clade3": clade2 & clade3,
}
for name, overlap in overlaps.items():
    if overlap:
        print(f"[WARN] Overlap detected in {name}: {len(overlap)} samples")

# 6) Build Unassigned = tree leaves not in any final clade
assigned = clade1a | clade1b | clade2 | clade3
unassigned = sorted([x for x in all_leaves if x.startswith("SRR") and x not in assigned])

with open(OUT_UNASSIGNED, "w") as f:
    f.write("\n".join(unassigned) + ("\n" if unassigned else ""))

# 7) Build metadata mapping
metadata_rows = []
for sample in sorted(all_leaves):
    if REF_NAME and sample == REF_NAME:
        lineage = "Reference"
    elif sample in clade1b:
        lineage = "Clade 1B"
    elif sample in clade1a:
        lineage = "Clade 1A"
    elif sample in clade2:
        lineage = "Clade 2"
    elif sample in clade3:
        lineage = "Clade 3"
    else:
        lineage = "Unassigned"

    metadata_rows.append({"Sample": sample, "Lineage": lineage})

df = pd.DataFrame(metadata_rows)

# 8) Save metadata
df.to_csv(OUT_METADATA, sep="\t", index=False)

# 9) Save counts summary
counts = df["Lineage"].value_counts(dropna=False)
with open(OUT_COUNTS, "w") as f:
    f.write("Lineage\tCount\n")
    for k, v in counts.items():
        f.write(f"{k}\t{v}\n")

# 10) Reporting
print(f"✅ Tree leaves extracted: {len(all_leaves)}")
print(f"✅ Clade 1 total: {len(clade1)}")
print(f"✅ Clade 1B total: {len(clade1b)}")
print(f"✅ Clade 1A total: {len(clade1a)}")
print(f"✅ Clade 2 total: {len(clade2)}")
print(f"✅ Clade 3 total: {len(clade3)}")
print(f"✅ Metadata written: {OUT_METADATA} ({len(df)} rows)")
print(f"✅ Unassigned written: {OUT_UNASSIGNED} ({len(unassigned)} SRR IDs)")
print(f"✅ Counts written: {OUT_COUNTS}")
print(counts)
