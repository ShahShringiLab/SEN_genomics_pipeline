#!/usr/bin/env python3
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]

checks = []

def require(path, tokens):
    text = (ROOT / path).read_text()
    for token in tokens:
        checks.append((f"{path}: {token}", token in text))

require("workflow/04_core_snp/13_snippy.sh", [
    "--mincov 10",
    "--minqual 100",
    "--mapqual 60",
    "--basequal 13",
    "--minfrac 0",
])

require("analysis/snps/19b_snps_byclade.py", [
    "FDR_ALPHA        = 0.05",
    "MIN_TOTAL_PRESENT = 10",
    "DEF_IN_CLADE     = 40.0",
    "DEF_OUTSIDE      = 1.0",
    'COVERAGE_BIN_ORDER = ["40-<60%", "60-<80%", "80-<90%", ">=90%"]',
    "total_k < M",
    'method="fdr_bh"',
])

require("workflow/06_assembly_pangenome/23_panaroo_byclade.py", [
    "ACCESSORY_MIN_FREQ = 0.01",
    "ACCESSORY_MAX_FREQ = 0.99",
    "MIN_TOTAL_PRESENT  = 10",
    "CLADEDEF_IN_PCT    = 90.0",
    "CLADEDEF_OUT_PCT   = 5.0",
    'method="fdr_bh"',
])

for path in [
    "analysis/clades/15_amr_byclade.py",
    "analysis/clades/16_plasmid_byclade.py",
    "analysis/clades/17_resfinder_byclade.py",
    "analysis/clades/18_vfdb_byclade.py",
]:
    require(path, [
        "MIN_COVERAGE = 80.0",
        "MIN_IDENTITY = 90.0",
        "MIN_TOTAL_PRESENT = 10",
        "MonteCarloMethod",
        "n_resamples=5000",
        "default_rng(12345)",
        'method="fdr_bh"',
        'method="holm"',
        "fisher_exact",
    ])

require("workflow/07_amr_vf_plasmid/25_abricate.sh", [
    'SEN_ABRICATE_MINCOV:-80',
    'SEN_ABRICATE_MINID:-90',
])

failures = [name for name, ok in checks if not ok]
for name, ok in checks:
    print(f"{'PASS' if ok else 'FAIL'}\t{name}")

if failures:
    print(f"\n[ERROR] {len(failures)} publication analysis-contract checks failed.", file=sys.stderr)
    sys.exit(1)

print(f"\n[PASS] {len(checks)} publication analysis-contract checks passed.")
