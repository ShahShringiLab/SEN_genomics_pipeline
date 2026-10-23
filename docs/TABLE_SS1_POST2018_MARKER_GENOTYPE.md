# Supplementary Table SS1 — Post-2018 QR/FQR marker genotype selection

This output supports the manuscript request for **each isolation year from 2018 through
2025**, stratified by QR/FQR marker genotype, including QS/FQS (no marker), with:

- **Available (N):** all matching genomes in the correctly parsed 33,831-record master.
- **Selected (pre-QC):** matching genomes among 3,434 archived isolate keys.
- **Analyzed (post-QC):** matching genomes among 3,306 final phylogeny Run IDs.

Run from repository root after activating the existing sen_stats environment:

```bash
bash scripts/run_master_part_a2.sh
```

The integrated master runner automatically exports:

- `results/part_a/master_33831/sampling_qc_table/Table_SS1_2018_2025_year_marker_genotype.xlsx`
- `results/part_a/master_33831/sampling_qc_table/Table_SS1_2018_2025_year_marker_genotype.csv`
- `results/part_a/master_33831/sampling_qc_table/Table_SS1_2018_2025_exact_genotype_signatures.csv`

The Excel file includes a concise manuscript table, an **exact marker-combination**
worksheet and a Definitions worksheet.

## Marker assignment and denominators

For the manuscript's concise table each genome is assigned to exactly one class:
`gyrA_D87Y`, `gyrA_D87N`, `gyrA_S83Y`, `gyrA_D87G`,
`qnrB19`, `Multiple major QR/FQR markers`, `Other QR/FQR marker(s)`,
or `QS/FQS (no marker)`. "Multiple major" applies to genomes carrying more than
one of the five specified major determinants. Additional markers on a genome that
carries exactly one major determinant are captured in its full exact genotype
signature but do not move it into "Multiple major".

This mutually exclusive grouping avoids double-counting and makes year totals
equal the numbers reported by the previously validated source x QR-status
annual sampling/QC table. Genotype calls come from the curated `AMR genotypes
core` field; QS/FQS here is **marker-negative**, not susceptibility confirmed
by phenotypic AST.

The source methods report only one representative from each outbreak cluster
per year and all 71 International isolates included in the selected cohort;
this export uses the **archived membership** rather than reapplying selection.

**Scope:** 2018 is included to match the 2018–2020 and 2021–2025
periods. These genotype-specific figures are NOT derivable from the year x source
x binary-QR summary alone: run the exporter with the entire parsed master.
