#!/usr/bin/env python3
# ============================================================
# HIGH-SPEED Prokka -> Panaroo pipeline (Threadripper-friendly)
#
# ✅ What this script does
#   1) Writes a reproducible environment setup script (mamba/conda)
#   2) Checks that the env + tools exist (prokka, panaroo, mafft)
#   3) Runs Prokka in parallel (many small jobs beats one big job)
#   4) Runs Panaroo using all cores
#
# ✅ Speed/robustness features
#   - Parallel Prokka with resume (skips genomes with existing .gff)
#   - Per-sample temp dirs to avoid collisions
#   - Suppresses noisy Prokka output (optional)
#   - Uses faster Prokka mode (--fast) for large N (3,000+)
#   - Panaroo clean-mode strict to reduce spurious cloud genes
#
# NOTE:
#   - This is optimized for speed, but the biggest bottleneck is I/O + Prokka
#   - For maximum throughput: use NVMe scratch, and avoid running on network FS
# ============================================================

import os
import glob
import subprocess
import multiprocessing
import sys
import shutil
from pathlib import Path
import textwrap

# -------------------------
# CONFIGURATION
# -------------------------
INPUT_DIR   = "./sistr_results_run/mini_assemblies"  # folder with .fasta/.fna/.fa
OUTPUT_DIR  = "./Panaroo_Run"
GFF_DIR     = os.path.join(OUTPUT_DIR, "gffs")

# Threadripper optimization:
CPU_CORES = multiprocessing.cpu_count()

# Prokka threads per process
PROKKA_CPUS_PER_JOB = 2

# Number of concurrent Prokka processes
# Heuristic: ~ half of cores / threads per job, but cap to avoid I/O thrash
# Example: 96 cores -> half is 48 jobs; with 2 threads/job => ~96 threads total
PROKKA_PROCESSES = max(1, int((CPU_CORES / 2) / max(PROKKA_CPUS_PER_JOB, 1)))

# Panaroo threads (typically all cores)
PANAROO_THREADS = CPU_CORES

# Prokka speed mode (good for 3000+ genomes; trades some sensitivity)
PROKKA_FAST = True

# Panaroo clean mode for large datasets
PANAROO_CLEAN_MODE = "strict"  # "strict" recommended for 3k+ genomes

# Create env helper file
ENV_NAME = "pangenome"
ENV_SCRIPT_PATH = os.path.join(OUTPUT_DIR, "00_create_env_pangenome.sh")

# -------------------------
# ENV CREATION SCRIPT WRITER
# -------------------------
def write_env_setup_script():
    """
    Writes a small bash script that creates a clean conda/mamba env with needed tools.
    This does NOT auto-run (safer). You run it once manually.
    """
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    script = textwrap.dedent(f"""\
    #!/usr/bin/env bash
    set -euo pipefail

    # ============================================================
    # Create high-performance pangenome env (Prokka + Panaroo + MAFFT)
    # Requires: mamba (recommended) or conda
    #
    # Usage:
    #   bash {Path(ENV_SCRIPT_PATH).name}
    #
    # Then:
    #   mamba activate {ENV_NAME}
    #   python {Path(__file__).name}
    # ============================================================

    if command -v mamba >/dev/null 2>&1; then
      PM="mamba"
    elif command -v conda >/dev/null 2>&1; then
      PM="conda"
    else
      echo "ERROR: conda or mamba not found. Install Miniforge/Mambaforge first."
      exit 1
    fi

    # Create env (idempotent)
    $PM create -y -n {ENV_NAME} -c conda-forge -c bioconda \\
      prokka \\
      panaroo \\
      mafft \\
      cd-hit \\
      blast \\
      hmmer \\
      parallel \\
      python=3.11 \\
      numpy pandas scipy networkx

    echo ""
    echo "✅ Env created: {ENV_NAME}"
    echo "Activate with: $PM activate {ENV_NAME}"
    """)

    with open(ENV_SCRIPT_PATH, "w") as f:
        f.write(script)
    os.chmod(ENV_SCRIPT_PATH, 0o755)

# -------------------------
# DEPENDENCY CHECK
# -------------------------
def check_dependencies():
    tools = ["prokka", "panaroo", "mafft"]
    missing = [t for t in tools if shutil.which(t) is None]

    if missing:
        print(f"❌ Missing tools in PATH: {', '.join(missing)}")
        print(f"➡️ I wrote an environment setup script you can run:")
        print(f"   {ENV_SCRIPT_PATH}")
        print("")
        print("Run:")
        print(f"  bash {ENV_SCRIPT_PATH}")
        print(f"  mamba activate {ENV_NAME}")
        print(f"  python {Path(__file__).name}")
        sys.exit(1)

    print(f"✅ Environment ready.")
    print(f"   Detected {CPU_CORES} CPU cores")
    print(f"   Prokka parallel jobs: {PROKKA_PROCESSES}  (each uses {PROKKA_CPUS_PER_JOB} threads)")
    print(f"   Panaroo threads: {PANAROO_THREADS}")

# -------------------------
# GENOME ANNOTATION WORKER
# -------------------------
def annotate_genome(file_path: str) -> str:
    """
    Runs Prokka for one genome, safely in parallel.
    - Resume-safe: if final GFF exists, skip
    - Writes to per-sample temp dir to prevent collisions
    """
    filename = os.path.basename(file_path)
    sample_name = os.path.splitext(filename)[0].replace("_trimmed", "")

    final_gff = os.path.join(GFF_DIR, f"{sample_name}.gff")
    if os.path.exists(final_gff) and os.path.getsize(final_gff) > 0:
        return f"SKIP\t{sample_name}\t(gff exists)"

    temp_root = os.path.join(OUTPUT_DIR, "temp_prokka")
    temp_out  = os.path.join(temp_root, sample_name)
    os.makedirs(temp_out, exist_ok=True)

    cmd = [
        "prokka",
        "--outdir", temp_out,
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

    cmd.append(file_path)

    try:
        # Send stdout/stderr to logs per-sample (faster to debug than DEVNULL)
        log_dir = os.path.join(OUTPUT_DIR, "logs_prokka")
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, f"{sample_name}.log")

        with open(log_file, "w") as lf:
            subprocess.run(cmd, stdout=lf, stderr=lf, check=True)

        src_gff = os.path.join(temp_out, f"{sample_name}.gff")
        if not os.path.exists(src_gff) or os.path.getsize(src_gff) == 0:
            return f"FAIL\t{sample_name}\t(no gff produced)"

        shutil.copy(src_gff, final_gff)

        # Cleanup per-sample temp folder (keep logs)
        shutil.rmtree(temp_out, ignore_errors=True)
        return f"OK\t{sample_name}"

    except subprocess.CalledProcessError:
        return f"FAIL\t{sample_name}\t(prokka error)"

# -------------------------
# MAIN
# -------------------------
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(GFF_DIR, exist_ok=True)
    os.makedirs(os.path.join(OUTPUT_DIR, "temp_prokka"), exist_ok=True)

    # Always write env helper (harmless)
    write_env_setup_script()

    # Check tools
    check_dependencies()

    # Gather FASTA files
    fasta_files = []
    for ext in ("*.fasta", "*.fna", "*.fa"):
        fasta_files.extend(glob.glob(os.path.join(INPUT_DIR, ext)))

    fasta_files = sorted(set(fasta_files))

    if not fasta_files:
        print(f"❌ No FASTA files found in: {INPUT_DIR}")
        sys.exit(1)

    print(f"📂 Found {len(fasta_files)} genomes.")

    # Prokka phase
    print("🚀 Starting Prokka annotation (parallel)...")
    # Use chunksize for speed (reduces IPC overhead)
    chunksize = max(1, len(fasta_files) // (PROKKA_PROCESSES * 10) if PROKKA_PROCESSES else 1)

    failed = 0
    with multiprocessing.Pool(processes=PROKKA_PROCESSES) as pool:
        for res in pool.imap_unordered(annotate_genome, fasta_files, chunksize=chunksize):
            # Print only failures (speed)
            if res.startswith("FAIL"):
                print(res)
                failed += 1

    print(f"✅ Prokka done. Failures: {failed}")

    # Remove temp root
    temp_root = os.path.join(OUTPUT_DIR, "temp_prokka")
    shutil.rmtree(temp_root, ignore_errors=True)

    # Panaroo phase
    print("🕸️ Starting Panaroo...")
    results_dir = os.path.join(OUTPUT_DIR, "results")
    os.makedirs(results_dir, exist_ok=True)

    # IMPORTANT: panaroo accepts glob if shell=True; otherwise expand list yourself
    gff_glob = os.path.join(GFF_DIR, "*.gff")
    cmd_panaroo = [
        "panaroo",
        "-i", gff_glob,
        "-o", results_dir,
        "--clean-mode", PANAROO_CLEAN_MODE,
        "--remove-invalid-genes",
        "-t", str(PANAROO_THREADS),
        "-a", "core",
        "--aligner", "mafft",
    ]

    # shell=True only because of wildcard expansion
    subprocess.run(" ".join(cmd_panaroo), shell=True, check=True)

    print(f"🎉 Finished! Results in: {results_dir}")
    print(f"🧪 Env helper script: {ENV_SCRIPT_PATH}")

if __name__ == "__main__":
    main()

