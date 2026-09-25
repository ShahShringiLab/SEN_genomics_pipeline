# Validation strategy

The publication-facing repository is validated in layers.

## Static portability and syntax

Run:

```bash
bash scripts/validate_repo.sh
```

This checks:

1. no workstation-specific paths or embedded Conda activation remain;
2. Bash syntax with `bash -n`;
3. Python syntax with `python3 -m py_compile`;
4. canonical metadata dimensions and final clade counts;
5. the current Git working-tree state.

## Tool provenance

Run the version capture script inside the relevant analysis environment(s):

```bash
bash scripts/capture_tool_versions.sh
```

The resulting `software_versions.tsv` is used to build reproducible environment
definitions without inventing package versions.

## Runtime validation

Static validation does not prove the biological workflow executes end-to-end.
After environments are resolved, run a small representative smoke cohort before
attempting the full 3,306-isolate reconstruction.
