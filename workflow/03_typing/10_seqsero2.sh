#!/bin/bash
# 1. Setup
source ~/miniforge3/etc/profile.d/conda.sh
conda activate seqsero2_env

IN="/home/samuelajulo/SENBio/Final/trimmed_fastq"
OUT="/home/samuelajulo/SENBio/Final/seqsero2_results"
mkdir -p "$OUT"

echo "=== STARTING SEQSERO2: $(date) ==="

# 2. Run with Parallel (Now that the tool exists!)
# -j 30 runs 30 samples at once
find "$IN" -name "*_trimmed_1.fastq.gz" | parallel -j 30 "
    R1={}
    R2=\${R1/_trimmed_1.fastq.gz/_trimmed_2.fastq.gz}
    SN=\$(basename \$R1 _trimmed_1.fastq.gz)
    
    SeqSero2_package.py -m k -t 2 -i \"\$R1\" \"\$R2\" -n \"\$SN\" -d \"$OUT/\$SN\" > /dev/null 2>&1
"

# 3. Merge Results
echo "=== MERGING RESULTS ==="
FIRST_FILE=$(find "$OUT" -name "SeqSero_result.tsv" | head -n 1)

if [ -n "$FIRST_FILE" ]; then
    head -n 1 "$FIRST_FILE" > "$OUT/SeqSero2_summary.tsv"
    find "$OUT" -name "SeqSero_result.tsv" -exec tail -n +2 {} + >> "$OUT/SeqSero2_summary.tsv"
    echo "SUCCESS: Results merged into $OUT/SeqSero2_summary.tsv"
else
    echo "ERROR: No results found. Did the installation work?"
fi

echo "=== FINISHED: $(date) ==="
