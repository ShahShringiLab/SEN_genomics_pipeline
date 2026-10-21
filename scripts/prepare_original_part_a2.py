#!/usr/bin/env python3
"""Prepare faithful original Part A2 Rmd from an uploaded notebook. No plot edits."""
from pathlib import Path
import hashlib
import sys
root=Path(__file__).resolve().parents[1]
original=Path(sys.argv[1]) if len(sys.argv)>1 else root/"provenance/SEN_Genomes2_uploaded_original_20261021.Rmd"
if not original.is_file(): raise SystemExit(f"Missing original Rmd: {original}")
text=original.read_text(encoding="utf-8-sig")
old='SEN <- read_csv(file_path, col_types = cols(.default = col_character()))'
assert text.count(old)==1, f"Expected exactly one original import, saw {text.count(old)}"
# Stage 0 parses and saves 33,831 x 51 before any derived columns. The notebook
# reads ONLY that snapshot; original analysis chunks remain verbatim.
new='''file_path <- file.path(Sys.getenv("SEN_ROOT", unset=getwd()),
                       "results/part_a/master_33831/00_parsed_raw/SEN_33831_PARSED_51COL.csv")
SEN <- read_csv(file_path, col_types = cols(.default = col_character()))
stopifnot(nrow(SEN)==33831L,ncol(SEN)==51L)'''
rendered=text.replace(old,new,1)
out=root/"analysis/part_a/original_rmd/SEN_Genomes2_PARSED_MASTER.Rmd"
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(rendered,encoding="utf-8")
# Strict invariant: only the input reader is changed.
assert rendered.replace(new,old,1)==text
sha=hashlib.sha256(original.read_bytes()).hexdigest()
(out.parent/"ORIGINAL_SHA256.txt").write_text(f"{sha}  {original.name}\n")
print(f"[PASS] Original notebook preserved; one import replaced: {out}")
print(f"[AUDIT] Original SHA-256: {sha}")
print("[NOTE] Notebook execution can depend on exploratory objects/interactive R chunks; validate outputs against the original.")
