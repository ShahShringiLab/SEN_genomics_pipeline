#!/bin/bash

# 1. Setup Paths
BASE_DIR="/home/samuelajulo/SENBio/Final"
ASSEMBLY_DIR="$BASE_DIR/shovill_assemblies"
ABRICATE_BASE_OUT="$BASE_DIR/abricate_results"

# ----------------------------------------------------
# CLEAN OUTPUT FOLDER BEFORE RUN (with safety guard)
# ----------------------------------------------------
# Guard 1: refuse if empty or root
if [[ -z "${ABRICATE_BASE_OUT:-}" || "$ABRICATE_BASE_OUT" == "/" ]]; then
    echo "[ERROR] ABRICATE_BASE_OUT is unsafe: '$ABRICATE_BASE_OUT' (refusing to rm -rf)"
    exit 1
fi
# Guard 2: refuse unless it lives under BASE_DIR
if [[ "$ABRICATE_BASE_OUT" != "$BASE_DIR/"* ]]; then
    echo "[ERROR] ABRICATE_BASE_OUT is not under BASE_DIR (refusing): '$ABRICATE_BASE_OUT'"
    exit 1
fi

if [ -d "$ABRICATE_BASE_OUT" ]; then
    echo "[INFO] Removing existing output folder: $ABRICATE_BASE_OUT"
    rm -rf "$ABRICATE_BASE_OUT"
fi
mkdir -p "$ABRICATE_BASE_OUT"
# ----------------------------------------------------

# Define the databases
DATABASES=("resfinder" "vfdb" "plasmidfinder")

for DB in "${DATABASES[@]}"; do
    echo "----------------------------------------------------"
    echo "Processing Database: $DB"
    DB_OUT_DIR="$ABRICATE_BASE_OUT/$DB"
    mkdir -p "$DB_OUT_DIR"
    
    # 2. Run Abricate
    # We use -j 16 to protect your IQ-TREE memory
    find "$ASSEMBLY_DIR" -name "contigs.fa" | parallel -j 16 "
        sample=\$(basename {//})
        output=\"$DB_OUT_DIR/\${sample}_\${DB}.tsv\"
        
        if [ ! -f \"\$output\" ]; then
            abricate --db $DB {} > \"\$output\"
        fi
    "

    # Wait for all background writes to finish
    sync

    # 3. Define Master Report INSIDE the subfolder
    MASTER_REPORT="$DB_OUT_DIR/master_${DB}_report.tsv"
    echo "Merging $DB results into $MASTER_REPORT..."

    # 4. Correct Merging Logic
    # Find the first .tsv file that IS NOT the master report itself
    FIRST_FILE=$(ls "$DB_OUT_DIR"/*.tsv 2>/dev/null | grep -v "master_" | head -n 1)
    
    if [ -z "$FIRST_FILE" ]; then
        echo "No individual result files found in $DB_OUT_DIR. Skipping merge."
        continue
    fi

    # Create Header with 'Sample' column
    echo -e "Sample\t$(head -n 1 "$FIRST_FILE")" > "$MASTER_REPORT"

    # Append data: Skip the header of each file and prepend Sample ID
    # Note: We grep -v to ensure we don't try to merge the master file into itself
    for f in "$DB_OUT_DIR"/*.tsv; do
        filename=$(basename "$f")
        if [[ "$filename" != master_* ]]; then
            sample=$(echo "$filename" | sed "s/_${DB}.tsv//")
            tail -n +2 "$f" | sed "s/^/${sample}\t/" >> "$MASTER_REPORT"
        fi
    done

    echo "Successfully created: $MASTER_REPORT"
done

echo "----------------------------------------------------"
echo "All tasks complete."
