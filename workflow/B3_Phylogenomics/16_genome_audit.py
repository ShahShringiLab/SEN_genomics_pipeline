#!/usr/bin/env python3

import os

import pandas as pd
import re
from pathlib import Path


# ============================================================
# CONFIGURATION
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[2]
BASE = Path(os.environ.get("SEN_ROOT", REPO_ROOT))

META_FILE = Path(os.environ.get("SEN_METADATA_FILE", BASE / "metadata" / "SEN_Genomes.csv"))
QC_FILE   = BASE / "three_step_QC_audit.csv"
TREE_FILE = BASE / "iqtree_final" / "SSLAB_FINAL.treefile"

OUTPUT_FILE = BASE / "SENGenomeAUDIT.csv"

MIN_COVERAGE = 30.0

# Expected historical cohort structure
EXPECTED_METADATA = 3434
EXPECTED_QC_PASS = 3309
EXPECTED_TREE_SRR = 3306
EXPECTED_TREE_TOTAL = 3307       # 3306 study SRRs + reference
EXPECTED_LOW_COV_RETAINED = 2
EXPECTED_IQTREE_EXCLUDED = 5


# ============================================================
# HELPERS
# ============================================================

def extract_srr(value):
    if pd.isna(value):
        return None

    m = re.search(
        r"SRR\d+",
        str(value),
        flags=re.IGNORECASE
    )

    return m.group(0).upper() if m else None


def clean_text(value):
    if pd.isna(value):
        return None

    x = str(value).strip()

    if x.lower() in {
        "",
        "nan",
        "na",
        "n/a",
        "none",
        "null"
    }:
        return None

    return x


def require_file(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Required file not found:\n{path}"
        )


# ============================================================
# VERIFY FILES
# ============================================================

for f in [
    META_FILE,
    QC_FILE,
    TREE_FILE
]:
    require_file(f)


print("=" * 88)
print("FINAL SALMONELLA ENTERITIDIS GENOME AUDIT")
print("=" * 88)

print(
    "\nHistorical cohort logic:\n"
    "  SEN_Genomes.csv\n"
    "       |\n"
    "       +-- Coverage / SISTR / SeqSero2 QC\n"
    "       |\n"
    "       +-- historical coverage exceptions where justified\n"
    "       |\n"
    "       +-- initial IQ-TREE composition/phylogenetic screening\n"
    "       |\n"
    "       +-- final IQ-TREE study cohort\n"
)


# ============================================================
# 1. ORIGINAL SEN_Genomes.csv
#
# Preserve every original column and row.
# Only AUDIT is appended.
# ============================================================

meta = pd.read_csv(
    META_FILE,
    dtype=str,
    keep_default_na=False,
    low_memory=False
)

meta.columns = [
    str(c).strip()
    for c in meta.columns
]

if "Run" not in meta.columns:
    raise RuntimeError(
        "SEN_Genomes.csv does not contain Run."
    )

original_columns = list(meta.columns)

meta["_AUDIT_SRR"] = (
    meta["Run"]
    .apply(extract_srr)
)

metadata_srrs = set(
    meta["_AUDIT_SRR"]
    .dropna()
)


print("-" * 88)
print("METADATA")
print("-" * 88)

print(
    f"SEN_Genomes.csv rows:                     "
    f"{len(meta):,}"
)

print(
    f"Unique metadata SRRs:                      "
    f"{len(metadata_srrs):,}"
)


# ============================================================
# 2. LOAD RETROSPECTIVE THREE-STEP QC
#
# This preserves EXACT 02g logic:
#
#   Coverage >=30X
#   SISTR Enteritidis
#   SeqSero2 Enteritidis
#
# We do NOT add new loci thresholds here.
# ============================================================

qc = pd.read_csv(
    QC_FILE,
    low_memory=False
)

qc.columns = [
    str(c).strip()
    for c in qc.columns
]

required_qc = [
    "Run",
    "FoldCoverage",
    "Coverage_Status",
    "SISTR_Call",
    "SISTR_Status",
    "SeqSero_Call",
    "SeqSero_Status",
    "Final_Status"
]

missing = [
    c for c in required_qc
    if c not in qc.columns
]

if missing:
    raise RuntimeError(
        "three_step_QC_audit.csv missing:\n"
        + "\n".join(missing)
    )


qc["Run"] = (
    qc["Run"]
    .apply(extract_srr)
)

qc["FoldCoverage"] = pd.to_numeric(
    qc["FoldCoverage"],
    errors="coerce"
)

qc = qc[
    qc["Run"].notna()
].copy()


# Ensure exactly one QC row per SRR
dups = qc[
    qc["Run"].duplicated(
        keep=False
    )
]

if not dups.empty:
    raise RuntimeError(
        "Duplicate SRRs found in three_step_QC_audit.csv:\n"
        + dups[
            ["Run", "Final_Status"]
        ].to_string(index=False)
    )


qc = qc.set_index(
    "Run",
    drop=False
)


three_step_pass = set(
    qc.loc[
        qc["Final_Status"] == "PASS_ALL_THREE",
        "Run"
    ]
)


print("\n" + "-" * 88)
print("THREE-STEP RETROSPECTIVE QC")
print("-" * 88)

print(
    f"QC records:                                "
    f"{len(qc):,}"
)

print(
    f"PASS_ALL_THREE:                            "
    f"{len(three_step_pass):,}"
)


# ============================================================
# 3. FINAL IQ-TREE
# ============================================================

tree_text = TREE_FILE.read_text()


tree_srrs = {
    x.upper()
    for x in re.findall(
        r"SRR\d+",
        tree_text,
        flags=re.IGNORECASE
    )
}


print("\n" + "-" * 88)
print("FINAL PHYLOGENY")
print("-" * 88)

print(
    f"SRR-labelled final-tree genomes:           "
    f"{len(tree_srrs):,}"
)

print(
    "Expected complete IQ-TREE:                 "
    "3,307 tips = 3,306 SRRs + 1 reference"
)


# ============================================================
# 4. IDENTIFY THE IMPORTANT HISTORICAL GROUPS
# ============================================================

both_serotype_pass = set(
    qc.loc[
        (qc["SISTR_Status"] == "PASS")
        &
        (qc["SeqSero_Status"] == "PASS"),
        "Run"
    ]
)


# ------------------------------------------------------------
# A. Retained final-tree genomes below retrospective 30X cutoff
# ------------------------------------------------------------

low_cov_retained = set(
    qc.loc[
        (qc["Run"].isin(tree_srrs))
        &
        (qc["SISTR_Status"] == "PASS")
        &
        (qc["SeqSero_Status"] == "PASS")
        &
        (qc["Coverage_Status"] == "FAIL"),
        "Run"
    ]
)


# ------------------------------------------------------------
# B. Passed all three retrospective filters, but absent final tree
#
# Historical downstream IQ-TREE exclusion.
# ------------------------------------------------------------

iqtree_excluded = set(
    qc.loc[
        (qc["Final_Status"] == "PASS_ALL_THREE")
        &
        (~qc["Run"].isin(tree_srrs)),
        "Run"
    ]
)


# ------------------------------------------------------------
# C. Ordinary final-tree genomes that passed all three
# ------------------------------------------------------------

standard_tree_pass = (
    tree_srrs
    &
    three_step_pass
)


print("\n" + "-" * 88)
print("SPECIAL COHORT GROUPS")
print("-" * 88)

print(
    f"Standard tree genomes passing all 3 QC:    "
    f"{len(standard_tree_pass):,}"
)

print(
    f"Below 30X but retained in final tree:       "
    f"{len(low_cov_retained):,}"
)

print(
    f"3-step PASS but excluded before final tree: "
    f"{len(iqtree_excluded):,}"
)


# ============================================================
# 5. FAILURE REASON FROM RETROSPECTIVE QC
# ============================================================

def qc_failure_reasons(row):

    reasons = []


    # --------------------------------------------------------
    # Coverage
    # --------------------------------------------------------

    cov_status = clean_text(
        row["Coverage_Status"]
    )

    coverage = row["FoldCoverage"]


    if cov_status == "FAIL":

        if pd.notna(coverage):

            reasons.append(
                f"Low Coverage "
                f"({coverage:.2f}X < {MIN_COVERAGE:.0f}X)"
            )

        else:

            reasons.append(
                "Coverage QC Failed"
            )


    elif cov_status == "MISSING":

        reasons.append(
            "Missing Coverage Data"
        )


    # --------------------------------------------------------
    # SISTR
    # --------------------------------------------------------

    sistr_status = clean_text(
        row["SISTR_Status"]
    )

    sistr_call = clean_text(
        row["SISTR_Call"]
    )


    if sistr_status == "FAIL":

        if sistr_call:

            reasons.append(
                f"SISTR Not Enteritidis "
                f"({sistr_call})"
            )

        else:

            reasons.append(
                "SISTR Not Enteritidis"
            )


    elif sistr_status == "MISSING":

        reasons.append(
            "Missing SISTR Data"
        )


    # --------------------------------------------------------
    # SeqSero2
    # --------------------------------------------------------

    seq_status = clean_text(
        row["SeqSero_Status"]
    )

    seq_call = clean_text(
        row["SeqSero_Call"]
    )


    if seq_status == "FAIL":

        if seq_call:

            reasons.append(
                f"SeqSero2 Not Enteritidis "
                f"({seq_call})"
            )

        else:

            reasons.append(
                "SeqSero2 Not Enteritidis"
            )


    elif seq_status == "MISSING":

        reasons.append(
            "Missing SeqSero2 Data"
        )


    return reasons


# ============================================================
# 6. BUILD FINAL AUDIT
# ============================================================

audit_map = {}


for run in sorted(metadata_srrs):


    # --------------------------------------------------------
    # No retrospective QC record
    # --------------------------------------------------------

    if run not in qc.index:

        audit_map[run] = (
            "FAIL: Missing from retrospective QC audit"
        )

        continue


    row = qc.loc[run]


    # ========================================================
    # CASE 1:
    # Genome is in the actual final IQ-TREE.
    #
    # Historical final inclusion takes precedence.
    # ========================================================

    if run in tree_srrs:


        # ----------------------------------------------------
        # Standard PASS
        # ----------------------------------------------------

        if row["Final_Status"] == "PASS_ALL_THREE":

            audit_map[run] = "PASS"

            continue


        # ----------------------------------------------------
        # Historical coverage exception:
        #
        # SISTR PASS
        # SeqSero2 PASS
        # Coverage just below retrospective 30X
        # Genome was nevertheless successfully retained in
        # the actual final phylogeny.
        # ----------------------------------------------------

        if (
            row["SISTR_Status"] == "PASS"
            and
            row["SeqSero_Status"] == "PASS"
            and
            row["Coverage_Status"] == "FAIL"
        ):

            cov = row["FoldCoverage"]

            audit_map[run] = (
                "PASS: Historical coverage exception; "
                "confirmed Enteritidis by SISTR and SeqSero2; "
                f"coverage {cov:.2f}X (<30X); "
                "retained in final IQ-TREE"
            )

            continue


        # ----------------------------------------------------
        # Unexpected:
        # final tree genome with some other retrospective QC
        # problem. Keep it visible.
        # ----------------------------------------------------

        reasons = qc_failure_reasons(
            row
        )

        audit_map[run] = (
            "PASS: Present in historical final IQ-TREE "
            "despite retrospective QC flag"
        )

        if reasons:

            audit_map[run] += (
                " ["
                + " | ".join(reasons)
                + "]"
            )

        continue


    # ========================================================
    # CASE 2:
    # Genome passed all retrospective three-step QC,
    # but is NOT in final tree.
    #
    # These are the historical IQ-TREE-screening exclusions.
    # ========================================================

    if row["Final_Status"] == "PASS_ALL_THREE":

        seq_call = clean_text(
            row["SeqSero_Call"]
        )


        # Explicitly preserve ambiguous antigenic call if one
        # exists, e.g. "Gallinarum or Enteritidis".
        ambiguous_seqsero = (
            seq_call is not None
            and
            re.search(
                r"\bor\b",
                seq_call,
                flags=re.IGNORECASE
            )
        )


        if ambiguous_seqsero:

            audit_map[run] = (
                "FAIL: Excluded during initial IQ-TREE "
                "composition/phylogenetic screening; "
                f"SeqSero2 ambiguous ({seq_call})"
            )

        else:

            audit_map[run] = (
                "FAIL: Excluded during initial IQ-TREE "
                "composition/phylogenetic screening"
            )

        continue


    # ========================================================
    # CASE 3:
    # Ordinary upstream QC failure.
    # ========================================================

    reasons = qc_failure_reasons(
        row
    )


    if reasons:

        audit_map[run] = (
            "FAIL: "
            + " | ".join(reasons)
        )

    else:

        audit_map[run] = (
            f"FAIL: Retrospective QC status "
            f"{row['Final_Status']}"
        )


# ============================================================
# 7. APPEND ONLY AUDIT TO ORIGINAL SEN_Genomes.csv
# ============================================================

def assign_audit(row):

    run = row["_AUDIT_SRR"]

    if run is None:

        return (
            "FAIL: Missing/Invalid Run accession"
        )

    return audit_map.get(
        run,
        "FAIL: Unable to resolve audit status"
    )


meta["AUDIT"] = (
    meta.apply(
        assign_audit,
        axis=1
    )
)


# EXACT original column order + AUDIT
meta = meta[
    original_columns
    + ["AUDIT"]
]


# ============================================================
# 8. WRITE MASTER AUDIT
# ============================================================

meta.to_csv(
    OUTPUT_FILE,
    index=False
)


# ============================================================
# 9. FINAL COUNTS
# ============================================================

audit_pass_srrs = {
    run
    for run, value in audit_map.items()
    if value.startswith("PASS")
}

audit_fail_srrs = {
    run
    for run, value in audit_map.items()
    if value.startswith("FAIL")
}


print("\n" + "=" * 88)
print("FINAL SEN GENOME AUDIT SUMMARY")
print("=" * 88)

print(
    f"Starting SEN genomes:                     "
    f"{len(metadata_srrs):,}"
)

print(
    f"Retrospective SISTR+SeqSero2 both PASS:    "
    f"{len(both_serotype_pass):,}"
)

print(
    f"Retrospective all-three QC PASS:           "
    f"{len(three_step_pass):,}"
)

print(
    f"Final study SRRs in IQ-TREE:               "
    f"{len(tree_srrs):,}"
)

print(
    f"Standard >=30X QC PASS in final tree:      "
    f"{len(standard_tree_pass):,}"
)

print(
    f"Historical <30X retained exceptions:       "
    f"{len(low_cov_retained):,}"
)

print(
    f"IQ-TREE screening exclusions:              "
    f"{len(iqtree_excluded):,}"
)

print(
    f"FINAL AUDIT PASS:                          "
    f"{len(audit_pass_srrs):,}"
)

print(
    f"FINAL AUDIT FAIL:                          "
    f"{len(audit_fail_srrs):,}"
)


# ============================================================
# 10. PRINT TWO HISTORICAL COVERAGE EXCEPTIONS
# ============================================================

print("\n" + "-" * 88)
print("HISTORICAL COVERAGE EXCEPTIONS RETAINED IN FINAL TREE")
print("-" * 88)

for run in sorted(
    low_cov_retained,
    key=lambda x: qc.loc[x, "FoldCoverage"]
):

    row = qc.loc[run]

    print(
        f"{run:15s} "
        f"{row['FoldCoverage']:7.2f}X  "
        f"SISTR={row['SISTR_Call']}  "
        f"SeqSero2={row['SeqSero_Call']}"
    )


# ============================================================
# 11. PRINT FIVE IQ-TREE EXCLUSIONS
# ============================================================

print("\n" + "-" * 88)
print("EXCLUDED DURING INITIAL IQ-TREE SCREENING")
print("-" * 88)

for run in sorted(
    iqtree_excluded,
    key=lambda x: qc.loc[x, "FoldCoverage"]
):

    row = qc.loc[run]

    print(
        f"{run:15s} "
        f"{row['FoldCoverage']:7.2f}X  "
        f"SISTR={row['SISTR_Call']}  "
        f"SeqSero2={row['SeqSero_Call']}"
    )


# ============================================================
# 12. HARD CONSISTENCY CHECKS
# ============================================================

print("\n" + "-" * 88)
print("COHORT CONSISTENCY CHECK")
print("-" * 88)

checks = []


def check(condition, success, failure):

    checks.append(condition)

    if condition:
        print(f"✓ {success}")

    else:
        print(f"✗ {failure}")


check(
    len(metadata_srrs) == EXPECTED_METADATA,
    "Starting cohort = 3,434",
    f"Starting cohort = {len(metadata_srrs):,}; expected 3,434"
)


check(
    len(three_step_pass) == EXPECTED_QC_PASS,
    "Retrospective all-three QC PASS = 3,309",
    f"All-three QC PASS = {len(three_step_pass):,}; expected 3,309"
)


check(
    len(tree_srrs) == EXPECTED_TREE_SRR,
    "Final phylogeny contains 3,306 SRR genomes",
    f"Final phylogeny contains {len(tree_srrs):,} SRRs; expected 3,306"
)


check(
    len(low_cov_retained) == EXPECTED_LOW_COV_RETAINED,
    "Exactly 2 historical <30X genomes were retained",
    f"Detected {len(low_cov_retained):,} retained <30X genomes; expected 2"
)


check(
    len(iqtree_excluded) == EXPECTED_IQTREE_EXCLUDED,
    "Exactly 5 >=30X three-step PASS genomes were excluded during IQ-TREE screening",
    f"Detected {len(iqtree_excluded):,} IQ-TREE exclusions; expected 5"
)


check(
    len(audit_pass_srrs) == EXPECTED_TREE_SRR,
    "FINAL AUDIT PASS = 3,306 study genomes",
    f"FINAL AUDIT PASS = {len(audit_pass_srrs):,}; expected 3,306"
)


check(
    audit_pass_srrs == tree_srrs,
    "AUDIT PASS set exactly matches final IQ-TREE SRR set",
    "AUDIT PASS set does NOT exactly match final IQ-TREE SRR set"
)


# ============================================================
# 13. AUDIT VALUE COUNTS
# ============================================================

print("\n" + "-" * 88)
print("AUDIT VALUE COUNTS")
print("-" * 88)

print(
    meta["AUDIT"]
    .value_counts()
    .to_string()
)


# ============================================================
# FINAL
# ============================================================

print("\n" + "=" * 88)

if all(checks):

    print(
        "✓ FINAL COHORT RECONSTRUCTION IS INTERNALLY CONSISTENT"
    )

    print(
        "\nFinal phylogeny:"
        "\n  3,306 study SRR genomes"
        "\n+     1 reference genome"
        "\n-------------------------"
        "\n  3,307 total IQ-TREE tips"
    )

else:

    print(
        "WARNING: One or more expected cohort counts "
        "did not reconcile."
    )


print(
    f"\nMaster audit written to:\n"
    f"  {OUTPUT_FILE}"
)

print("=" * 88)

