#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/../lib/common.sh"

require_cmd SeqSero2_package.py
require_cmd parallel

IN="${SEN_TRIMMED_READS}"
OUT="${SEN_SEQSERO2_OUT:-$SEN_ROOT/seqsero2_results}"
JOBS="${SEN_SEQSERO2_JOBS:-30}"

mkdir -p "$OUT"

echo "=== STARTING SEQSERO2: $(date) ==="

find "$IN" -name "*_trimmed_1.fastq.gz" | parallel -j "$JOBS" "
    R1={}
    R2=\${R1/_trimmed_1.fastq.gz/_trimmed_2.fastq.gz}
    SN=\$(basename \$R1 _trimmed_1.fastq.gz)

    SeqSero2_package.py -m k -t 2 -i \"\$R1\" \"\$R2\" -n \"\$SN\" -d \"$OUT/\$SN\" > /dev/null 2>&1
"

echo "=== MERGING RESULTS ==="
FIRST_FILE=$(find "$OUT" -name "SeqSero_result.tsv" | head -n 1)

if [[ -n "$FIRST_FILE" ]]; then
    head -n 1 "$FIRST_FILE" > "$OUT/SeqSero2_summary.tsv"
    find "$OUT" -name "SeqSero_result.tsv" -exec tail -n +2 {} + >> "$OUT/SeqSero2_summary.tsv"
    echo "SUCCESS: Results merged into $OUT/SeqSero2_summary.tsv"
else
    echo "ERROR: No results found."
    exit 1
fi

echo "=== FINISHED: $(date) ==="
