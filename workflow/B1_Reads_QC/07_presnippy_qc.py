import pandas as pd
import numpy as np

# --- CONFIGURATION ---
INPUT_CSV = "pilot_full_db_results.csv"
OUTPUT_CSV = "final_qc_annotated_3k.csv"
WHITELIST_FILE = "snippy_ready_srrs.txt"

# Salmonella Enteritidis Parameters
GENOME_SIZE = 4800000  # 4.8 Mbp
MIN_PURITY = 95.0      # % Salmonella
MIN_DEPTH = 30.0       # 30x coverage
MIN_READ_LEN = 75.0    # 75bp (Standard for PE150 quality)

def run_qc_pipeline():
    print(f"🚀 Loading results from {INPUT_CSV}...")
    try:
        df = pd.read_csv(INPUT_CSV)
    except FileNotFoundError:
        print("❌ Error: pilot_full_db_results.csv not found yet.")
        return

    # --- 1. METRIC CALCULATIONS ---
    
    # Total Bases * 2 (because SeqKit only scanned R1, and we have paired ends)
    df['Estimated_Depth'] = (df['Total_Bases'] * 2) / GENOME_SIZE
    
    # Average length of a single read after trimming
    df['Mean_Read_Length'] = df['Total_Bases'] / df['Total_Reads']
    
    # Purity is already in your 'Salmonella_Percent' column

    # --- 2. MULTI-METRIC FILTERING ---
    
    def apply_filters(row):
        reasons = []
        if row['Salmonella_Percent'] < MIN_PURITY:
            reasons.append("Low_Purity")
        if row['Estimated_Depth'] < MIN_DEPTH:
            reasons.append("Low_Depth")
        if row['Mean_Read_Length'] < MIN_READ_LEN:
            reasons.append("Short_Reads")
        
        return "PASS" if not reasons else "|".join(reasons)

    df['QC_Status'] = df.apply(apply_filters, axis=1)

    # --- 3. SAVE DATA & WHITELIST ---
    
    df.to_csv(OUTPUT_CSV, index=False)
    
    # Create the Whitelist (Only PASSing SRRs)
    whitelist = df[df['QC_Status'] == "PASS"]['Sample']
    whitelist.to_csv(WHITELIST_FILE, index=False, header=False)

    # --- 4. SUMMARY REPORT ---
    pass_count = len(whitelist)
    total_count = len(df)
    pass_rate = (pass_count / total_count) * 100

    print("\n" + "="*40)
    print("      GENOMIC QC DASHBOARD")
    print("="*40)
    print(f"Total Genomes Analyzed:   {total_count}")
    print(f"Genomes Passing QC:       {pass_count} ({pass_rate:.2f}%)")
    print("-" * 40)
    print("Failure Breakdown:")
    print(df['QC_Status'].value_counts().drop('PASS', errors='ignore'))
    print("-" * 40)
    print(f"✅ Whitelist saved to: {WHITELIST_FILE}")
    print(f"📊 Full report saved to: {OUTPUT_CSV}")
    print("="*40)

if __name__ == "__main__":
    run_qc_pipeline()
