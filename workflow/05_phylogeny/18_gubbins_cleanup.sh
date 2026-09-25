# 1. Create the exclusion list
cat <<EOF > failed_samples.txt
SRR15168826
SRR34428946
SRR3634462
SRR5459684
EOF

# 2. Filter the FASTA (using awk for safety)
# This keeps only sequences whose headers ARE NOT in the failed_samples file
awk 'BEGIN{while((getline < "failed_samples.txt") > 0) f[">"$1]=1} /^>/ {skip=(f[$1])} !skip' \
  gubbin/senbio_res.filtered_polymorphic_sites.fasta > gubbin/clean_final_alignment.fasta

# 3. Verify the count (should be 3307)
echo "Original count: $(grep -c ">" gubbin/senbio_res.filtered_polymorphic_sites.fasta)"
echo "Cleaned count: $(grep -c ">" gubbin/clean_final_alignment.fasta)"

# 4. Kickstart IQ-TREE (Update your 11_iqtree.sh to use clean_final_alignment.fasta first!)
bash 11_iqtree.sh
