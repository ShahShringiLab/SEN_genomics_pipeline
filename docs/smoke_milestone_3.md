# Smoke milestone 3 — assembly, annotation and pangenome

Status: **PASS**

Representative cohort: four *Salmonella Enteritidis* isolates.

Validated stages:

1. Shovill draft assembly from Kraken-cleaned paired reads
2. Prokka annotation configured for *Salmonella enterica*
3. Panaroo 1.6.0 strict-mode pangenome construction
4. MAFFT core-genome alignment

Observed smoke result:

- Shovill assemblies: 4/4
- Prokka GFFs: 4/4
- Panaroo gene presence/absence matrix: produced successfully
- Panaroo core alignment: produced successfully
- Panaroo graph processing completed on 4,456 elements
- 4,285 of 4,344 core genes were retained after Panaroo filtering

This smoke milestone validates runtime and stage contracts only. It does not
constitute accessory-gene enrichment analysis for the full study.
