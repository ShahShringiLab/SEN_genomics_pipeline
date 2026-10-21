#!/usr/bin/env python3
"""Generate R execution script from original Rmd chunks 4-21, without editing code."""
from pathlib import Path
import re,hashlib
root=Path(__file__).resolve().parents[1]
source=root/"analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
if not source.is_file(): raise SystemExit("First run bash scripts/run_master_part_a2.sh")
src=source.read_text(encoding="utf-8")
parts=re.findall(r'(?ms)^```\{r[^\n}]*\}\s*\n(.*?)^```\s*$',src)
assert len(parts)==26,f"Expected exactly 26 original chunks, observed {len(parts)}"
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
for index in range(4,22):
    preamble.extend([f'\n# BEGIN ORIGINAL CHUNK {index}',parts[index-1],
       guard%(index,index),f'# END ORIGINAL CHUNK {index}'])
rfile.write_text("\n".join(preamble)+"\n")
print(f"[READY] Original A2 chunks 4-21, with row integrity checks: {rfile}")
