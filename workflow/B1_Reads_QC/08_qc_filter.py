#!/usr/bin/env python3

import os

import pandas as pd
import re
import sys
from pathlib import Path


# ============================================================
# CONFIGURATION
# ============================================================

REPO_ROOT = Path(__file__).resolve().parents[2]
BASE = Path(os.environ.get("SEN_ROOT", REPO_ROOT))

META_FILE = Path(os.environ.get("SEN_METADATA_FILE", BASE / "metadata" / "SEN_Genomes.csv"))

COV_FILE = (
    BASE /
    "trimmed_fastq_qc" /
    "all_samples_coverage.csv"
)

SISTR_FILE = (
    BASE /
    "sistr_results_run" /
    "sistr_master_summary.csv"
)

SEQSERO_FILE = (
    BASE /
    "seqsero2_results" /
    "SeqSero2_summary.tsv"
)

MIN_COVERAGE = 30.0


# ============================================================
# OUTPUTS
# ============================================================

OUT_AUDIT = (
    BASE /
    "three_step_QC_audit.csv"
)

OUT_PASS_ALL = (
    BASE /
    "all_three_pass_srrs.txt"
)

OUT_FAIL_ALL = (
    BASE /
    "all_three_fail_srrs.txt"
)

OUT_COVERAGE_FAIL = (
    BASE /
    "coverage_fail_srrs.txt"
)

OUT_SISTR_FAIL = (
    BASE /
    "sistr_fail_srrs.txt"
)

OUT_SEQSERO_FAIL = (
    BASE /
    "seqsero_fail_srrs.txt"
)

OUT_COV_SISTR_FAIL = (
    BASE /
    "coverage_and_sistr_fail_srrs.txt"
)

OUT_COV_SEQSERO_FAIL = (
    BASE /
    "coverage_and_seqsero_fail_srrs.txt"
)

OUT_SISTR_SEQSERO_FAIL = (
    BASE /
    "sistr_and_seqsero_fail_srrs.txt"
)

OUT_ANY_FAIL = (
    BASE /
    "any_qc_fail_srrs.txt"
)

OUT_MISSING = (
    BASE /
    "qc_missing_data_srrs.txt"
)


# ============================================================
# HELPERS
# ============================================================

def extract_srr(value):
    """
    Extract SRR accession from arbitrary text/path.

    Examples:
        SRR123456
        SRR123456.fasta
        /path/to/SRR123456.fasta
        SRR123456_trimmed_1.fastq.gz
    """

    if pd.isna(value):
        return None

    match = re.search(
        r"SRR\d+",
        str(value),
        flags=re.IGNORECASE
    )

    if match:
        return match.group(0).upper()

    return None


def clean_columns(df):
    df.columns = [
        str(c).strip()
        for c in df.columns
    ]
    return df


def save_srr_list(series, outfile):
    """
    Save sorted unique SRR list without header.
    """

    values = (
        series
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
    )

    values.to_csv(
        outfile,
        index=False,
        header=False
    )

    return len(values)


def is_enteritidis(value):
    """
    Require the word Enteritidis in the relevant serotype field.
    """

    if pd.isna(value):
        return False

    return bool(
        re.search(
            r"\bEnteritidis\b",
            str(value),
            flags=re.IGNORECASE
        )
    )


# ============================================================
# INPUT CHECK
# ============================================================

print("=" * 76)
print("SEN THREE-STEP QC FILTER")
print("=" * 76)

required_files = [
    META_FILE,
    COV_FILE,
    SISTR_FILE,
    SEQSERO_FILE,
]

missing_files = [
    str(x)
    for x in required_files
    if not x.exists()
]

if missing_files:

    print("\nERROR: Required input file(s) missing:")

    for x in missing_files:
        print(f"  {x}")

    sys.exit(1)


# ============================================================
# STEP 0
# SEN_Genomes.csv IS THE DENOMINATOR
# ============================================================

print("\n" + "=" * 76)
print("STEP 0: SEN GENOME UNIVERSE")
print("=" * 76)

meta = pd.read_csv(
    META_FILE,
    low_memory=False
)

meta = clean_columns(meta)

if "Run" not in meta.columns:

    raise RuntimeError(
        "SEN_Genomes.csv does not contain a 'Run' column."
    )


meta["Run"] = (
    meta["Run"]
    .apply(extract_srr)
)


bad_meta = meta[
    meta["Run"].isna()
]

if not bad_meta.empty:

    print(
        f"WARNING: {len(bad_meta):,} metadata row(s) "
        "do not contain a valid SRR accession."
    )


universe = (
    meta[
        ["Run"]
    ]
    .dropna()
    .drop_duplicates()
    .sort_values("Run")
    .reset_index(drop=True)
)


print(
    f"Unique SRRs in SEN_Genomes.csv: "
    f"{len(universe):,}"
)


if universe.empty:

    raise RuntimeError(
        "SEN_Genomes.csv produced zero SRRs."
    )


# ============================================================
# STEP 1
# COVERAGE >= 30X
# ============================================================

print("\n" + "=" * 76)
print("STEP 1: COVERAGE QC")
print("=" * 76)

coverage = pd.read_csv(
    COV_FILE,
    low_memory=False
)

coverage = clean_columns(
    coverage
)


if "FoldCoverage" not in coverage.columns:

    raise RuntimeError(
        "Coverage file does not contain 'FoldCoverage'.\n"
        f"Available columns:\n{list(coverage.columns)}"
    )


# ------------------------------------------------------------
# Find likely sample identifier column
# ------------------------------------------------------------

coverage_id_candidates = [
    "Sample",
    "sample",
    "sample_id",
    "Run",
    "run",
    "File",
    "file",
    "Filename",
    "filename",
]


coverage_id_col = None

for candidate in coverage_id_candidates:

    if candidate in coverage.columns:
        coverage_id_col = candidate
        break


if coverage_id_col is None:

    coverage_id_col = coverage.columns[0]

    print(
        "WARNING: Could not find a standard coverage ID column."
    )

    print(
        f"Using first column instead: {coverage_id_col}"
    )


print(
    f"Coverage sample-ID column: "
    f"{coverage_id_col}"
)


coverage["Run"] = (
    coverage[
        coverage_id_col
    ]
    .apply(extract_srr)
)


coverage["FoldCoverage"] = pd.to_numeric(
    coverage["FoldCoverage"],
    errors="coerce"
)


# ------------------------------------------------------------
# Check duplicated coverage SRRs
# ------------------------------------------------------------

coverage_valid = (
    coverage[
        coverage["Run"].notna()
    ]
    .copy()
)


duplicate_cov = (
    coverage_valid[
        coverage_valid["Run"].duplicated(
            keep=False
        )
    ]
)


if not duplicate_cov.empty:

    print(
        f"WARNING: {duplicate_cov['Run'].nunique():,} "
        "SRR(s) have multiple coverage rows."
    )

    # Conservative / deterministic:
    # retain maximum observed FoldCoverage per SRR
    coverage_small = (
        coverage_valid
        .groupby(
            "Run",
            as_index=False
        )[
            "FoldCoverage"
        ]
        .max()
    )

else:

    coverage_small = (
        coverage_valid[
            [
                "Run",
                "FoldCoverage"
            ]
        ]
        .drop_duplicates(
            subset=["Run"]
        )
    )


universe = universe.merge(
    coverage_small,
    on="Run",
    how="left",
    validate="1:1"
)


universe["Coverage_Status"] = "MISSING"


universe.loc[
    universe["FoldCoverage"] >= MIN_COVERAGE,
    "Coverage_Status"
] = "PASS"


universe.loc[
    universe["FoldCoverage"] < MIN_COVERAGE,
    "Coverage_Status"
] = "FAIL"


print(
    universe[
        "Coverage_Status"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


# ============================================================
# STEP 2
# SISTR ENTERITIDIS
# ============================================================

print("\n" + "=" * 76)
print("STEP 2: SISTR ENTERITIDIS")
print("=" * 76)

sistr = pd.read_csv(
    SISTR_FILE,
    low_memory=False
)

sistr = clean_columns(
    sistr
)


if "genome" not in sistr.columns:

    raise RuntimeError(
        "SISTR file does not contain 'genome'.\n"
        f"Available columns:\n{list(sistr.columns)}"
    )


sistr["Run"] = (
    sistr["genome"]
    .apply(extract_srr)
)


# ------------------------------------------------------------
# Prefer cgMLST serovar when available
# ------------------------------------------------------------

if "serovar_cgmlst" in sistr.columns:

    SISTR_CALL_COL = "serovar_cgmlst"

elif "serovar" in sistr.columns:

    SISTR_CALL_COL = "serovar"

else:

    raise RuntimeError(
        "SISTR file contains neither "
        "'serovar_cgmlst' nor 'serovar'."
    )


print(
    f"SISTR classification column: "
    f"{SISTR_CALL_COL}"
)


sistr[
    "SISTR_Call"
] = sistr[
    SISTR_CALL_COL
].astype(str)


sistr[
    "SISTR_Enteritidis"
] = sistr[
    SISTR_CALL_COL
].apply(
    is_enteritidis
)


# ------------------------------------------------------------
# Duplicate/conflict audit
# ------------------------------------------------------------

sistr_valid = (
    sistr[
        sistr["Run"].notna()
    ]
    .copy()
)


duplicate_sistr = (
    sistr_valid[
        sistr_valid["Run"].duplicated(
            keep=False
        )
    ]
)


if not duplicate_sistr.empty:

    conflicts = (
        duplicate_sistr
        .groupby(
            "Run"
        )[
            "SISTR_Enteritidis"
        ]
        .nunique()
    )

    conflicts = conflicts[
        conflicts > 1
    ]

    if not conflicts.empty:

        print(
            "\nERROR: Conflicting SISTR calls detected:"
        )

        print(
            duplicate_sistr[
                duplicate_sistr[
                    "Run"
                ].isin(
                    conflicts.index
                )
            ][
                [
                    "Run",
                    SISTR_CALL_COL,
                    "genome"
                ]
            ]
            .sort_values("Run")
            .to_string(index=False)
        )

        raise RuntimeError(
            "Resolve conflicting SISTR calls first."
        )


sistr_small = (
    sistr_valid[
        [
            "Run",
            "SISTR_Call",
            "SISTR_Enteritidis"
        ]
    ]
    .drop_duplicates(
        subset=["Run"]
    )
)


universe = universe.merge(
    sistr_small,
    on="Run",
    how="left",
    validate="1:1"
)


universe["SISTR_Status"] = "MISSING"


universe.loc[
    universe["SISTR_Enteritidis"] == True,
    "SISTR_Status"
] = "PASS"


universe.loc[
    universe["SISTR_Enteritidis"] == False,
    "SISTR_Status"
] = "FAIL"


print(
    universe[
        "SISTR_Status"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


# ============================================================
# STEP 3
# SEQSERO2 ENTERITIDIS
# ============================================================

print("\n" + "=" * 76)
print("STEP 3: SEQSERO2 ENTERITIDIS")
print("=" * 76)

seqsero = pd.read_csv(
    SEQSERO_FILE,
    sep="\t",
    low_memory=False
)

seqsero = clean_columns(
    seqsero
)


required_seqsero_columns = [
    "Sample name",
    "Predicted serotype",
]


missing_seqsero_columns = [
    c
    for c in required_seqsero_columns
    if c not in seqsero.columns
]


if missing_seqsero_columns:

    raise RuntimeError(
        "SeqSero2 file missing required columns:\n"
        + "\n".join(
            missing_seqsero_columns
        )
        + "\n\nAvailable columns:\n"
        + str(
            list(seqsero.columns)
        )
    )


# ------------------------------------------------------------
# CLEAN SEQSERO2 SUMMARY
#
# The combined SeqSero2 summary may contain shell/cat separator
# rows such as:
#
# ==> /path/to/SRR9984402/SeqSero_result.tsv <==
#
# These are NOT biological results and must be removed before
# duplicate/conflict assessment.
# ------------------------------------------------------------

seqsero["Sample name"] = (
    seqsero["Sample name"]
    .astype(str)
    .str.strip()
)

seqsero["Predicted serotype"] = (
    seqsero["Predicted serotype"]
    .astype(str)
    .str.strip()
)

# Flag synthetic file-header/separator rows
separator_mask = (
    seqsero["Sample name"].str.contains(
        r"==>|<==",
        regex=True,
        na=False
    )
)

print(
    f"SeqSero2 synthetic separator/header rows removed: "
    f"{separator_mask.sum():,}"
)

seqsero = (
    seqsero[
        ~separator_mask
    ]
    .copy()
)

# Extract the SRR only from genuine SeqSero2 result rows
seqsero["Run"] = (
    seqsero[
        "Sample name"
    ]
    .apply(extract_srr)
)

# Treat empty/nan-like serotype values as missing, not FAIL
seqsero["SeqSero_Call"] = (
    seqsero[
        "Predicted serotype"
    ]
    .replace(
        {
            "nan": pd.NA,
            "NaN": pd.NA,
            "NA": pd.NA,
            "N/A": pd.NA,
            "": pd.NA,
        }
    )
)

seqsero["SeqSero_Enteritidis"] = pd.NA

valid_seqsero_call = (
    seqsero["SeqSero_Call"].notna()
)

seqsero.loc[
    valid_seqsero_call,
    "SeqSero_Enteritidis"
] = (
    seqsero.loc[
        valid_seqsero_call,
        "SeqSero_Call"
    ]
    .apply(
        is_enteritidis
    )
)


# ------------------------------------------------------------
# Duplicate/conflict audit
# ------------------------------------------------------------

seqsero_valid = (
    seqsero[
        seqsero["Run"].notna()
    ]
    .copy()
)


duplicate_seqsero = (
    seqsero_valid[
        seqsero_valid["Run"].duplicated(
            keep=False
        )
    ]
)


if not duplicate_seqsero.empty:

    conflicts = (
        duplicate_seqsero
        .groupby(
            "Run"
        )[
            "SeqSero_Enteritidis"
        ]
        .nunique()
    )

    conflicts = conflicts[
        conflicts > 1
    ]

    if not conflicts.empty:

        print(
            "\nERROR: Conflicting SeqSero2 calls detected:"
        )

        print(
            duplicate_seqsero[
                duplicate_seqsero[
                    "Run"
                ].isin(
                    conflicts.index
                )
            ][
                [
                    "Run",
                    "Sample name",
                    "Predicted serotype"
                ]
            ]
            .sort_values("Run")
            .to_string(index=False)
        )

        raise RuntimeError(
            "Resolve conflicting SeqSero2 calls first."
        )


seqsero_small = (
    seqsero_valid[
        [
            "Run",
            "SeqSero_Call",
            "SeqSero_Enteritidis"
        ]
    ]
    .drop_duplicates(
        subset=["Run"]
    )
)


universe = universe.merge(
    seqsero_small,
    on="Run",
    how="left",
    validate="1:1"
)


universe[
    "SeqSero_Status"
] = "MISSING"


universe.loc[
    universe[
        "SeqSero_Enteritidis"
    ] == True,
    "SeqSero_Status"
] = "PASS"


universe.loc[
    universe[
        "SeqSero_Enteritidis"
    ].eq(False),
    "SeqSero_Status"
] = "FAIL"


print(
    universe[
        "SeqSero_Status"
    ]
    .value_counts(
        dropna=False
    )
    .to_string()
)


# ============================================================
# FINAL THREE-STEP STATUS
# ============================================================

print("\n" + "=" * 76)
print("FINAL CLASSIFICATION")
print("=" * 76)


def classify_final(row):

    cov = row[
        "Coverage_Status"
    ]

    sis = row[
        "SISTR_Status"
    ]

    seq = row[
        "SeqSero_Status"
    ]


    # --------------------------------------------------------
    # Missing information handled separately
    # --------------------------------------------------------

    missing = []

    if cov == "MISSING":
        missing.append("COVERAGE")

    if sis == "MISSING":
        missing.append("SISTR")

    if seq == "MISSING":
        missing.append("SEQSERO")


    if missing:

        return (
            "MISSING_"
            + "_".join(missing)
        )


    # --------------------------------------------------------
    # All pass
    # --------------------------------------------------------

    if (
        cov == "PASS"
        and sis == "PASS"
        and seq == "PASS"
    ):

        return "PASS_ALL_THREE"


    # --------------------------------------------------------
    # Determine actual failing filters
    # --------------------------------------------------------

    failed = []

    if cov == "FAIL":
        failed.append("COVERAGE")

    if sis == "FAIL":
        failed.append("SISTR")

    if seq == "FAIL":
        failed.append("SEQSERO")


    if len(failed) == 3:
        return "FAIL_ALL_THREE"

    if len(failed) == 2:
        return (
            "FAIL_"
            + "_".join(failed)
        )

    if len(failed) == 1:
        return (
            "FAIL_"
            + failed[0]
            + "_ONLY"
        )

    return "CHECK"


universe[
    "Final_Status"
] = universe.apply(
    classify_final,
    axis=1
)


print(
    universe[
        "Final_Status"
    ]
    .value_counts()
    .to_string()
)


# ============================================================
# BOOLEAN CONVENIENCE COLUMNS
# ============================================================

universe[
    "Coverage_PASS"
] = (
    universe[
        "Coverage_Status"
    ] == "PASS"
)


universe[
    "SISTR_PASS"
] = (
    universe[
        "SISTR_Status"
    ] == "PASS"
)


universe[
    "SeqSero_PASS"
] = (
    universe[
        "SeqSero_Status"
    ] == "PASS"
)


universe[
    "PASS_ALL"
] = (
    universe[
        "Final_Status"
    ] == "PASS_ALL_THREE"
)


# ============================================================
# SAVE COMPLETE AUDIT
# ============================================================

universe.to_csv(
    OUT_AUDIT,
    index=False
)


# ============================================================
# EXTRACT INDIVIDUAL FAILURES
# ============================================================

coverage_fail = universe.loc[
    universe[
        "Coverage_Status"
    ] == "FAIL",
    "Run"
]


sistr_fail = universe.loc[
    universe[
        "SISTR_Status"
    ] == "FAIL",
    "Run"
]


seqsero_fail = universe.loc[
    universe[
        "SeqSero_Status"
    ] == "FAIL",
    "Run"
]


# ============================================================
# PAIRWISE FAILURE LISTS
#
# These mean that BOTH named criteria failed,
# regardless of the third criterion.
# ============================================================

coverage_sistr_fail = universe.loc[
    (
        universe[
            "Coverage_Status"
        ] == "FAIL"
    )
    &
    (
        universe[
            "SISTR_Status"
        ] == "FAIL"
    ),
    "Run"
]


coverage_seqsero_fail = universe.loc[
    (
        universe[
            "Coverage_Status"
        ] == "FAIL"
    )
    &
    (
        universe[
            "SeqSero_Status"
        ] == "FAIL"
    ),
    "Run"
]


sistr_seqsero_fail = universe.loc[
    (
        universe[
            "SISTR_Status"
        ] == "FAIL"
    )
    &
    (
        universe[
            "SeqSero_Status"
        ] == "FAIL"
    ),
    "Run"
]


# ============================================================
# ALL THREE PASS / FAIL
# ============================================================

pass_all = universe.loc[
    universe[
        "Final_Status"
    ] == "PASS_ALL_THREE",
    "Run"
]


fail_all = universe.loc[
    universe[
        "Final_Status"
    ] == "FAIL_ALL_THREE",
    "Run"
]


# ============================================================
# ANY ACTUAL QC FAILURE
# ============================================================

any_fail = universe.loc[
    (
        universe[
            "Coverage_Status"
        ] == "FAIL"
    )
    |
    (
        universe[
            "SISTR_Status"
        ] == "FAIL"
    )
    |
    (
        universe[
            "SeqSero_Status"
        ] == "FAIL"
    ),
    "Run"
]


# ============================================================
# MISSING DATA
# ============================================================

missing_data = universe.loc[
    (
        universe[
            "Coverage_Status"
        ] == "MISSING"
    )
    |
    (
        universe[
            "SISTR_Status"
        ] == "MISSING"
    )
    |
    (
        universe[
            "SeqSero_Status"
        ] == "MISSING"
    ),
    "Run"
]


# ============================================================
# WRITE LISTS
# ============================================================

n_cov_fail = save_srr_list(
    coverage_fail,
    OUT_COVERAGE_FAIL
)

n_sistr_fail = save_srr_list(
    sistr_fail,
    OUT_SISTR_FAIL
)

n_seqsero_fail = save_srr_list(
    seqsero_fail,
    OUT_SEQSERO_FAIL
)

n_cov_sistr = save_srr_list(
    coverage_sistr_fail,
    OUT_COV_SISTR_FAIL
)

n_cov_seq = save_srr_list(
    coverage_seqsero_fail,
    OUT_COV_SEQSERO_FAIL
)

n_sistr_seq = save_srr_list(
    sistr_seqsero_fail,
    OUT_SISTR_SEQSERO_FAIL
)

n_pass_all = save_srr_list(
    pass_all,
    OUT_PASS_ALL
)

n_fail_all = save_srr_list(
    fail_all,
    OUT_FAIL_ALL
)

n_any_fail = save_srr_list(
    any_fail,
    OUT_ANY_FAIL
)

n_missing = save_srr_list(
    missing_data,
    OUT_MISSING
)


# ============================================================
# FINAL REPORT
# ============================================================

print("\n" + "=" * 76)
print("THREE-STEP QC SUMMARY")
print("=" * 76)

print(
    f"Total SEN genomes:                 "
    f"{len(universe):,}"
)

print(
    f"PASS all three:                    "
    f"{n_pass_all:,}"
)

print(
    f"Coverage failures:                 "
    f"{n_cov_fail:,}"
)

print(
    f"SISTR failures:                    "
    f"{n_sistr_fail:,}"
)

print(
    f"SeqSero2 failures:                 "
    f"{n_seqsero_fail:,}"
)

print(
    f"Coverage + SISTR failures:         "
    f"{n_cov_sistr:,}"
)

print(
    f"Coverage + SeqSero2 failures:      "
    f"{n_cov_seq:,}"
)

print(
    f"SISTR + SeqSero2 failures:         "
    f"{n_sistr_seq:,}"
)

print(
    f"FAIL all three:                    "
    f"{n_fail_all:,}"
)

print(
    f"Any true QC failure:               "
    f"{n_any_fail:,}"
)

print(
    f"Missing >=1 QC result:             "
    f"{n_missing:,}"
)


# ============================================================
# PRINT ALL-THREE FAILURES
# ============================================================

print("\n" + "-" * 76)
print("SRRs FAILING ALL THREE FILTERS")
print("-" * 76)

if n_fail_all == 0:

    print("None")

else:

    for run in (
        fail_all
        .drop_duplicates()
        .sort_values()
    ):

        print(run)


# ============================================================
# FILE REPORT
# ============================================================

print("\n" + "=" * 76)
print("FILES WRITTEN")
print("=" * 76)

print(
    f"Full audit:\n"
    f"  {OUT_AUDIT}"
)

print(
    f"\nPASS all three:\n"
    f"  {OUT_PASS_ALL}"
)

print(
    f"\nFAIL all three:\n"
    f"  {OUT_FAIL_ALL}"
)

print(
    f"\nCoverage failures:\n"
    f"  {OUT_COVERAGE_FAIL}"
)

print(
    f"\nSISTR failures:\n"
    f"  {OUT_SISTR_FAIL}"
)

print(
    f"\nSeqSero2 failures:\n"
    f"  {OUT_SEQSERO_FAIL}"
)

print(
    f"\nCoverage + SISTR failures:\n"
    f"  {OUT_COV_SISTR_FAIL}"
)

print(
    f"\nCoverage + SeqSero2 failures:\n"
    f"  {OUT_COV_SEQSERO_FAIL}"
)

print(
    f"\nSISTR + SeqSero2 failures:\n"
    f"  {OUT_SISTR_SEQSERO_FAIL}"
)

print(
    f"\nAny QC failure:\n"
    f"  {OUT_ANY_FAIL}"
)

print(
    f"\nMissing QC data:\n"
    f"  {OUT_MISSING}"
)

print("\n" + "=" * 76)
print("DONE")
print("=" * 76)

