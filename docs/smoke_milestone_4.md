# Smoke milestone 4 — AMR, virulence and plasmid screening

Status: **PASS**

Representative cohort: four *Salmonella Enteritidis* assemblies.

Validated stages:

1. AMRFinderPlus 4.2.5 with `--plus --organism Salmonella`
2. frozen AMRFinderPlus database release `2026-08-07.1`
3. ABRicate 1.2.0 ResFinder screening
4. ABRicate 1.2.0 VFDB screening
5. ABRicate 1.2.0 PlasmidFinder screening
6. manuscript presence thresholds encoded at screening time:
   coverage >=80%, identity >=90%

Observed database inventory:

| Database | Sequences | Type | Date |
|---|---:|---|---|
| ResFinder | 3,206 | nucl | 2025-Dec-5 |
| VFDB | 4,592 | nucl | 2025-Dec-5 |
| PlasmidFinder | 488 | nucl | 2025-Dec-5 |

Observed smoke result:

- AMRFinderPlus reports: 4/4
- ResFinder reports: 4/4
- VFDB reports: 4/4
- PlasmidFinder reports: 4/4
- end-to-end smoke pipeline: PASS

The repository now pins the AMRFinder database by exact dated release rather
than downloading whichever database is current at runtime. The exact database
source path is reconstructed from NCBI's published AMRFinderPlus database
layout and validated after local indexing.
