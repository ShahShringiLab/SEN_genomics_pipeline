#!/usr/bin/env python3
"""Generate R execution from master-population Rmd chunks, excluding cohort selection."""
from pathlib import Path
import re,hashlib
root=Path(__file__).resolve().parents[1]
source=root/"analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
if not source.is_file(): raise SystemExit("First run bash scripts/run_master_part_a2.sh")
src=source.read_text(encoding="utf-8")
parts=re.findall(r'(?ms)^```\{r[^\n}]*\}\s*\n(.*?)^```\s*$',src)
assert len(parts) in (26, 30),f"Expected 26 or 30 original chunks, observed {len(parts)}"
selection_ix=[i for i,chunk in enumerate(parts,1) if
    "FINAL 3500-GENOME SELECTION SCRIPT" in chunk]
assert len(selection_ix)==1,f"Cannot uniquely locate cohort selection chunk: {selection_ix}"
last_master_chunk=selection_ix[0]-1
assert 21 <= last_master_chunk <= 22, (
    f"Unexpected master-analysis chunk boundary: {last_master_chunk}"
)

out=root/"results/part_a/master_33831/original_execution"
out.mkdir(parents=True,exist_ok=True)
rfile=out/"run_original_A2_master.R"
preamble=[
'options(warn=1)',
'library(dplyr); library(tidyr); library(readr); library(stringr)',
'library(ggplot2); library(lubridate); library(purrr); library(patchwork)',
'library(scales); library(forcats)',
f'cat("[INFO] Original notebook source sha256: {hashlib.sha256(src.encode()).hexdigest()}\\n")',
]
guard='''if (!exists("SEN") || nrow(SEN) != 33831L ||
    dplyr::n_distinct(SEN$Run) != 33831L)
  stop("Master 33,831 population integrity failed after chunk %d")
cat("[PASS] Original A2 chunk %d; 33,831 genomes remain\\n")
'''
for index in range(4,last_master_chunk+1):
    preamble.extend([f'\n# BEGIN ORIGINAL CHUNK {index}',parts[index-1],
       guard%(index,index),f'# END ORIGINAL CHUNK {index}'])
    if "DOTPLOT BY YEAR" in parts[index-1]:
        preamble.extend([
            'if (!exists("den_year")) stop("Temporal dot plot denominator missing")',
            'if (sum(den_year$N) != 33508L) stop(paste0("Temporal dot plot N expected 33,508, got ",sum(den_year$N)))',
            'cat("[PASS] Temporal dot plot: 33,508 dated genomes, including all years <=2000\\n")',
        ])
rfile.write_text("\n".join(preamble)+"\n")
print(f"[READY] Original A2 master chunks 4-{last_master_chunk} (from {len(parts)} Rmd chunks), with 33,831-row and 33,508-temporal integrity checks: {rfile}")
