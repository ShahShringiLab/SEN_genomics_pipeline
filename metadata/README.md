# Metadata

Canonical inputs for reconstruction of the analysis cohort live here.

## Files

- `SEN_Genomes.csv` — original 3,434-genome selected cohort metadata. This
  file must be copied unchanged from `SENBio-Final@78abf1e`.
- `final_clade_metadata.tsv` — canonical final tree-derived clade assignment,
  imported from `ITOL/Clade_metadata.txt`.
- `clade_counts.tsv` — checksum of final clade membership.

The final clade metadata contains 3,307 rows: 3,306 analyzed isolates and the
P125109 reference.

Expected final counts:

| Lineage | Count |
|---|---:|
| Clade 3 | 1357 |
| Clade 2 | 894 |
| Clade 1B | 759 |
| Clade 1A | 169 |
| Unassigned | 127 |
| Reference | 1 |

The similarly named historical `Clade_metadata_updated.txt` is not the
complete canonical input because it excludes Unassigned isolates and the
reference.
