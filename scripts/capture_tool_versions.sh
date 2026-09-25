#!/usr/bin/env bash
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/software_versions.tsv}"

printf "tool\tversion_output\n" > "$OUT"

capture() {
  local name="$1"
  shift
  local out
  out="$("$@" 2>&1 | head -n 1 || true)"
  out="${out//$'\t'/ }"
  out="${out//$'\n'/ }"
  printf "%s\t%s\n" "$name" "$out" >> "$OUT"
}

capture prefetch prefetch --version
capture fasterq-dump fasterq-dump --version
capture fastp fastp --version
capture fastqc fastqc --version
capture multiqc multiqc --version
capture kraken2 kraken2 --version
capture seqkit seqkit version
capture SeqSero2 SeqSero2_package.py --version
capture skesa skesa --version
capture sistr sistr --version
capture mlst mlst --version
capture snippy snippy --version
capture snippy-core snippy-core --version
capture snpEff snpEff -version
capture gubbins run_gubbins.py --version
capture veryfasttree veryfasttree -version
capture iqtree2 iqtree2 --version
capture shovill shovill --version
capture spades spades.py --version
capture prokka prokka --version
capture panaroo panaroo --version
capture amrfinder amrfinder --version
capture abricate abricate --version
capture parallel parallel --version
capture python3 python3 --version
capture R R --version

echo "Wrote: $OUT"
