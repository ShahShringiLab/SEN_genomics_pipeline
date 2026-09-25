#!/usr/bin/env python3

import glob
import multiprocessing
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SEN_ROOT = Path(os.environ.get("SEN_ROOT", REPO_ROOT))

INPUT_DIR = Path(os.environ.get("SEN_FINAL_CONTIGS_DIR", SEN_ROOT / "Final_Contigs_Only"))
OUTPUT_DIR = Path(os.environ.get("SEN_PROKKA_OUT", SEN_ROOT / "prokka_annotations"))
LOG_DIR = Path(os.environ.get("SEN_PROKKA_LOG_DIR", SEN_ROOT / "logs" / "prokka"))

CPU_CORES = multiprocessing.cpu_count()
PROKKA_CPUS_PER_JOB = int(os.environ.get("SEN_PROKKA_CPUS_PER_JOB", "2"))
PROKKA_PROCESSES = int(
    os.environ.get(
        "SEN_PROKKA_JOBS",
        str(max(1, int((CPU_CORES / 2) / max(PROKKA_CPUS_PER_JOB, 1))))
    )
)
PROKKA_FAST = os.environ.get("SEN_PROKKA_FAST", "1") == "1"

def require_cmd(cmd: str):
    if shutil.which(cmd) is None:
        raise SystemExit(f"[ERROR] Required command not found in PATH: {cmd}")

def annotate_genome(file_path: str) -> str:
    fp = Path(file_path)
    sample_name = fp.stem.replace("_trimmed", "")
    final_gff = OUTPUT_DIR / f"{sample_name}.gff"

    if final_gff.exists() and final_gff.stat().st_size > 0:
        return f"SKIP\t{sample_name}\t(gff exists)"

    temp_out = OUTPUT_DIR / "_tmp" / sample_name
    temp_out.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    cmd = [
        "prokka",
        "--outdir", str(temp_out),
        "--prefix", sample_name,
        "--locustag", sample_name,
        "--genus", "Salmonella",
        "--species", "enterica",
        "--strain", sample_name,
        "--cpus", str(PROKKA_CPUS_PER_JOB),
        "--force",
    ]
    if PROKKA_FAST:
        cmd.append("--fast")
    cmd.append(str(fp))

    log_file = LOG_DIR / f"{sample_name}.log"

    try:
        with log_file.open("w") as lf:
            subprocess.run(cmd, stdout=lf, stderr=lf, check=True)

        src_gff = temp_out / f"{sample_name}.gff"
        if not src_gff.exists() or src_gff.stat().st_size == 0:
            return f"FAIL\t{sample_name}\t(no gff produced)"

        shutil.copy2(src_gff, final_gff)
        shutil.rmtree(temp_out, ignore_errors=True)
        return f"OK\t{sample_name}"
    except subprocess.CalledProcessError:
        return f"FAIL\t{sample_name}\t(prokka error)"

def main():
    require_cmd("prokka")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "_tmp").mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    fasta_files = []
    for ext in ("*.fasta", "*.fna", "*.fa"):
        fasta_files.extend(glob.glob(str(INPUT_DIR / ext)))
    fasta_files = sorted(set(fasta_files))

    if not fasta_files:
        raise SystemExit(f"[ERROR] No FASTA files found in {INPUT_DIR}")

    print(f"[INFO] Found {len(fasta_files)} Shovill assemblies")
    print(f"[INFO] Prokka jobs={PROKKA_PROCESSES}, cpus/job={PROKKA_CPUS_PER_JOB}")
    print(f"[INFO] --fast={'enabled' if PROKKA_FAST else 'disabled'}")

    failed = 0
    chunksize = max(1, len(fasta_files) // (PROKKA_PROCESSES * 10))

    with multiprocessing.Pool(processes=PROKKA_PROCESSES) as pool:
        for result in pool.imap_unordered(annotate_genome, fasta_files, chunksize=chunksize):
            if result.startswith("FAIL"):
                print(result)
                failed += 1

    shutil.rmtree(OUTPUT_DIR / "_tmp", ignore_errors=True)

    if failed:
        raise SystemExit(f"[ERROR] Prokka failures: {failed}")

    print(f"[INFO] Prokka complete: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
