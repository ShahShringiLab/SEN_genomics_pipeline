# Pipeline provenance

This publication-facing repository is derived from the verified analysis
snapshot:

- Source repository: `ajulojays/SENBio-Final`
- Frozen source commit: `78abf1e0efb3ed2eb91d7de72fc928a3b90e64f5`
- Short SHA: `78abf1e`
- Snapshot status at freeze: no modified tracked source files; remaining
  workstation-only files were untracked generated/runtime outputs.

The source snapshot provides the analysis provenance record. Files in this
repository are organized, parameterized, documented, and tested components of
the complete publication workflow rather than a verbatim workstation mirror.

## Provenance policy

1. Preserve scientific logic from the frozen source snapshot.
2. Preserve canonical metadata separately from generated outputs.
3. Remove hard-coded workstation paths.
4. Move machine-specific database locations into local configuration.
5. Keep generated results out of version control.
6. Record release-level implementation changes explicitly in version control.
