import os
import glob
import pandas as pd
import re
from pathlib import Path

# 1. SETUP
REPO_ROOT = Path(__file__).resolve().parents[2]
BASE = Path(os.environ.get("SEN_ROOT", REPO_ROOT))
OUT_DIR = Path(os.environ.get("SEN_SNIPPY_OUT", BASE / "Snippy_output"))
FINAL_REPORT = OUT_DIR / "SENBIO_RECOVERED_REPORT.csv"

def parse_snpeff_vcf(vcf_path, sample_id):
    """
    Extracts POS, REF, ALT, GENE, and AA_CHANGE from SnpEff VCF.
    Handles the ANN field format: 
    Allele | Annotation | Impact | GeneName | GeneID | FeatureType | FeatureID | ... | HGVS.p
    """
    rows = []
    if not os.path.exists(vcf_path):
        return rows

    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            cols = line.strip().split('\t')
            if len(cols) < 8:
                continue
                
            pos, ref, alt, info = cols[1], cols[3], cols[4], cols[7]
            
            # Look for the ANN field created by SnpEff
            match = re.search(r'ANN=([^;]+)', info)
            if match:
                # SnpEff can list multiple effects separated by commas; we take the first (most significant)
                first_ann = match.group(1).split(',')[0]
                ann_data = first_ann.split('|')
                
                # EFFECT is index 1
                effect = ann_data[1] if len(ann_data) > 1 else "unknown"
                
                # GENE NAME logic: Use Gene Name (index 3), fallback to Locus Tag/ID (index 4)
                gene_name = ann_data[3] if len(ann_data) > 3 else ""
                gene_id = ann_data[4] if len(ann_data) > 4 else ""
                gene = gene_name if gene_name else (gene_id if gene_id else "intergenic")
                
                # AA_CHANGE (HGVS.p) is index 10
                aa_change = ann_data[10] if len(ann_data) > 10 else "."
                if not aa_change: aa_change = "."
            else:
                gene, effect, aa_change = "unannotated", "unknown", "."
            
            rows.append({
                'Sample': sample_id,
                'POS': int(pos),
                'REF': ref,
                'ALT': alt,
                'GENE': gene,
                'EFFECT': effect,
                'AA_CHANGE': aa_change
            })
    return rows

# 2. PROCESS ALL SAMPLES
all_variants = []
# Find all annotated VCFs in the subdirectories
vcf_files = glob.glob(os.path.join(OUT_DIR, "*", "snps.annotated.vcf"))

if not vcf_files:
    print(f"❌ No 'snps.annotated.vcf' files found in {OUT_DIR}")
else:
    print(f"🚀 Found {len(vcf_files)} annotated VCFs. Starting extraction...")

    for vcf in vcf_files:
        # Get the folder name as the Sample ID
        sample_name = os.path.basename(os.path.dirname(vcf))
        all_variants.extend(parse_snpeff_vcf(vcf, sample_name))

    # 3. CREATE FINAL DATAFRAME
    if all_variants:
        df = pd.DataFrame(all_variants)

        # Quick summary stats
        total_muts = len(df)
        # Filter for gyrA (case-insensitive check)
        gyra_df = df[df['GENE'].str.contains('gyrA', case=False, na=False)]
        gyra_count = len(gyra_df)
        
        print("-" * 30)
        print(f"Total mutations extracted: {total_muts}")
        print(f"✅ gyrA mutations found: {gyra_count}")
        
        # Show top 5 gyrA mutations if any exist
        if gyra_count > 0:
            print("Preview of gyrA hits:")
            print(gyra_df[['Sample', 'POS', 'AA_CHANGE']].head())

        # Save to CSV
        df.to_csv(FINAL_REPORT, index=False)
        print("-" * 30)
        print(f"🎉 Success! Master report saved to: {FINAL_REPORT}")
    else:
        print("Empty VCFs or no variants found.")
