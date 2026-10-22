#!/usr/bin/env python3
"""Fast preflight: verify immutable 30-chunk original SEN Rmd and extracted workflow.

No 33,831-genome R calculation is executed here. It is safe before a full rerun.
"""
import hashlib
from pathlib import Path
import re
import subprocess
import sys

root = Path(__file__).resolve().parents[1]
original = root / "provenance/SEN_Genomes2_uploaded_original_20261021.Rmd"
if not original.is_file():
    raise SystemExit("[FAIL] Complete original user-uploaded Rmd missing")
s = original.read_text(encoding="utf-8-sig")
sha = hashlib.sha256(s.encode()).hexdigest()
EXPECTED_SHA = "3f4c1451ebbe8fe603efd514a087dc60f841d15157096169ef70749ca92d51b2"
if sha != EXPECTED_SHA:
    raise SystemExit(
        f"[FAIL] Original Rmd source changed: observed sha256={sha}, expected={EXPECTED_SHA}")
chunks = re.findall(r'(?ms)^\x60\x60\x60\{r[^\n}]*\}\s*\n(.*?)^\x60\x60\x60\s*$', s)
if len(chunks) != 30:
    raise SystemExit(f"[FAIL] Original Rmd expected 30 R chunks, found {len(chunks)}")
if "FINAL 3500-GENOME SELECTION SCRIPT" not in chunks[22]:
    raise SystemExit("[FAIL] Historical selection must occur at chunk 23")
if "SEN_TableS1_final_source_decomposition.csv" not in chunks[16]:
    raise SystemExit("[FAIL] Full-master source-refinement chunk 17 missing")
if "DOTPLOT BY YEAR" not in chunks[18]:
    raise SystemExit("[FAIL] Full-master temporal chart chunk 19 missing")

subprocess.run([sys.executable, str(root / "scripts/prepare_original_part_a2.py"),
                str(original)], cwd=root, check=True)
subprocess.run([sys.executable, str(root / "scripts/build_original_a2_master_execution.py")],
               cwd=root, check=True)
working = (root / "analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd").read_text()
master = (root / "results/part_a/master_33831/original_execution/run_original_A2_master.R").read_text()
combined = (root / "results/part_a/master_33831/original_execution/run_original_A2_all.R").read_text()
checks = [
    ('30 source chunks preserved in working Rmd',
     len(re.findall(r'(?ms)^\x60\x60\x60\{r[^\n}]*\}\s*\n(.*?)^\x60\x60\x60\s*$',working)) == 30),
    ('original source preserved', original.read_text(encoding="utf-8-sig") == s),
    ('horizontal marker chart reaches 250,000, with endpoint tick omitted',
     'limits = c(0, 250000)' in working and
     'breaks = seq(0, 200000, by = 50000)' in working),
    ('historical year filters removed from working full-master',
     'Collection_year >= 1990,' not in "\n".join(chunks2 for chunks2 in
        re.findall(r'(?ms)^\x60\x60\x60\{r[^\n}]*\}\s*\n(.*?)^\x60\x60\x60\s*$', working)[3:22])),
    ('original full-master chunks source table', 'SEN_TableS1_final_source_decomposition.csv' in master),
    ('main full-master input check', '33831L' in master),
    ('all full-master Rmd-derived fields are exported',
     'SEN_33831_original_Rmd_derived.csv' in master),
    ('full-master stage never uses selected QC input', 'Finalupdated_SEN_Genomes1_QC_merged.csv' not in master),
    ('separate selected original stage uses QC input', 'Finalupdated_SEN_Genomes1_QC_merged.csv' in combined),
    ('selected plotting inputs contain QR_status and corrected QR genotype checks',
     combined.count('QR_status=dplyr::if_else(QR_genotype==1L')==2
     and 'SEN$QR_genotype <- as.integer(SEN$FQR_genotype)' in combined),
    ('selected resume script was generated from original master snapshot',
     (root / "results/part_a/master_33831/original_execution/run_original_A2_selected_resume.R").is_file()),
    ('historical selection generation skipped', 'BEGIN USER ORIGINAL CHUNK 23' not in combined),
]
for label, ok in checks:
    if not ok:
        raise SystemExit(f"[FAIL] Preflight: {label}")
    print(f"[PASS] Preflight: {label}")
print("[PASS] Original-Rmd static preflight; ready for actual 33,831-genome R execution")
