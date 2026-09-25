#!/usr/bin/env bash
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$ROOT/software_versions.tsv}"

printf "tool\tstatus\texecutable\tversion_output\n" > "$OUT"

probe() {
  local name="$1"
  shift

  if ! command -v "$1" >/dev/null 2>&1; then
    printf "%s\tNOT_FOUND\t\t\n" "$name" >> "$OUT"
    return
  fi

  local exe out
  exe="$(command -v "$1")"
  out="$("$@" 2>&1 | sed '/^[[:space:]]*$/d' | head -n 1 || true)"
  out="${out//$'\t'/ }"
  out="${out//$'\n'/ }"
  printf "%s\tFOUND\t%s\t%s\n" "$name" "$exe" "$out" >> "$OUT"
}

probe prefetch prefetch --version
probe fasterq-dump fasterq-dump --version
probe fastp fastp --version
probe fastqc fastqc --version
probe multiqc multiqc --version
probe kraken2 kraken2 --version
probe seqkit seqkit version
probe SeqSero2 SeqSero2_package.py --version
probe skesa skesa --version
probe sistr sistr --version
probe mlst mlst --version
probe snippy snippy --version
probe snippy-core snippy-core --version
probe snpEff snpEff -version
probe gubbins run_gubbins.py --version
probe veryfasttree veryfasttree -version
probe iqtree2 iqtree2 --version
probe shovill shovill --version
probe spades spades.py --version
probe prokka prokka --version
probe panaroo panaroo --version
probe amrfinder amrfinder --version
probe abricate abricate --version
probe parallel parallel --version
probe python3 python3 --version
probe R R --version

echo "Wrote: $OUT"
