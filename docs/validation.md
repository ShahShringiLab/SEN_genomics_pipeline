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

The resulting `software_versions.tsv` records the resolved runtime software
versions for reproducibility and release verification.

## Runtime validation

Part A runtime validation is performed with:

```bash
bash scripts/run_part_a.sh
```

The runner validates committed inputs and the Part A -> Part B accession handoff
before execution, runs all five Part A R stages, and then requires
`scripts/validate_part_a_outputs.py` to pass before emitting `SEN PART A:
PASS`.

Part B runtime validation uses the representative four-isolate smoke cohort to
verify end-to-end stage handoffs before manuscript-scale execution.
