# Smoke milestone 2 — core SNPs, recombination filtering and phylogeny

Status: **PASS**

Representative cohort: four *Salmonella Enteritidis* isolates plus the P125109
reference chromosome (NC_011294.1).

Validated stages:

1. automatic P125109 reference bootstrap and validation
2. Snippy on all four isolates
3. snippy-core alignment
4. Gubbins recombination filtering
5. historical final-tree exclusion step
6. IQ-TREE 2.2.6 with ASC-aware model search, 1000 UFBoot2 replicates and BNNI

Observed smoke result:

- Snippy isolate outputs: 4/4
- core alignment taxa: 5/5
- Gubbins filtered-alignment taxa: 5/5
- final cleaned-alignment taxa: 5/5
- IQ-TREE maximum-likelihood tree: produced successfully
- IQ-TREE UFBoot samples: 1000

The four-isolate smoke test validates executable workflow contracts only. It is
not intended to reproduce the topology or clade structure of the full
3,306-isolate study.
