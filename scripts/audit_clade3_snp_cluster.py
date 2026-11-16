#!/usr/bin/env python3
"""Independent, standard-library check of Clade 3/NCBI SNP-cluster manuscript numbers.

Inputs are frozen SEN source files in the GitHub repository:
- data/part_a/genotype/US-SEN_12192025.csv (33,831 master isolates)
- data/part_a/selection/Finalupdated_SEN_Genomes1.csv (3,434 pre-QC selected)
- metadata/final_clade_metadata.tsv (3,306 post-QC isolates)

Uses the exact in-memory CSV syntax normalization documented by
workflow/A2_Genotype/lib/common.R. No source files are modified.
"""
import csv
import io
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGET = "PDS000065758"
ISOLATE_ID = re.compile(r"^SRR\d+$")
SUFFIX = re.compile(r"\.\d+$")


def canonical_cluster(x):
    return SUFFIX.sub("", (x or "").strip())


def d87y_positive(amr):
    return any(
        token.strip().split("=", 1)[0].strip().lower() == "gyra_d87y"
        for token in (amr or "").split(",")
    )


def read_csv(path, normalize=False):
    raw = path.read_text(encoding="utf-8-sig")
    if normalize:
        assert "__SEN_LITERAL_COMMA__" not in raw
        n_esc = raw.count("\\,")
        n_quote = raw.count(',""Exudate, Left Heel"","clinical",')
        raw = raw.replace("\\,", "__SEN_LITERAL_COMMA__")
        raw = raw.replace(',""Exudate, Left Heel"","clinical",',
                          ',"Exudate, Left Heel","clinical",')
        print(f"[PARSER] normalized escaped commas={n_esc}; quoted Exudate fields={n_quote}")
    reader = csv.DictReader(io.StringIO(raw, newline=""))
    headers = reader.fieldnames
    if not headers:
        raise ValueError(f"No header in {path}")
    rows = []
    for lineno, row in enumerate(reader, start=2):
        if None in row or any(value is None for value in row.values()):
            raise ValueError(f"Malformed CSV row {lineno} in {path}: {row}")
        if normalize:
            row = {key: value.replace("__SEN_LITERAL_COMMA__", ",") for key, value in row.items()}
        rows.append(row)
    print(f"[FILE] {path.relative_to(ROOT)}: rows={len(rows)} columns={len(headers)}")
    return rows, headers


def main():
    full, master_columns = read_csv(ROOT / "data/part_a/genotype/US-SEN_12192025.csv", normalize=True)
    assert len(full) == 33831, f"Expected 33831, got {len(full)}"
    assert len(master_columns) == 51, f"Expected 51 columns, got {len(master_columns)}"
    assert len({r["Run"] for r in full}) == 33831, "Duplicate Run in full master"
    assert "SNP cluster" in master_columns
    assert "AMR genotypes core" in master_columns

    selected, selected_columns = read_csv(ROOT / "data/part_a/selection/Finalupdated_SEN_Genomes1.csv")
    assert len(selected) == 3434
    assert "SNP cluster" in selected_columns
    assert "Run" in selected_columns
    clade_rows, _ = read_csv(ROOT / "metadata/final_clade_metadata.tsv") if False else (None, None)
    with (ROOT / "metadata/final_clade_metadata.tsv").open(newline="") as f:
        clade_reader = csv.DictReader(f, delimiter="\t")
        clade_rows = list(clade_reader)
    clade_map = {r["Sample"]: r["Lineage"] for r in clade_rows if r["Sample"] != "Reference"}
    assert len(clade_map) == 3306
    assert all(ISOLATE_ID.fullmatch(x) for x in clade_map)

    selected_study = [r for r in selected if r["Run"] in clade_map]
    assert len(selected_study) == 3306
    clade3 = [r for r in selected_study if clade_map[r["Run"]] == "Clade 3"]
    target_c3 = [r for r in clade3 if canonical_cluster(r["SNP cluster"]) == TARGET]
    target_else = [r for r in selected_study
                   if canonical_cluster(r["SNP cluster"]) == TARGET and clade_map[r["Run"]] != "Clade 3"]
    assert len(clade3) == 1357
    assert len(target_c3) == 1343
    assert len(target_else) == 0
    print(f"[CLADE 3] matching={len(target_c3)}/{len(clade3)} ({100*len(target_c3)/len(clade3):.4f}%)")
    print(f"[CLADE 3] target cluster in other selected clades={len(target_else)}")

    positives = [r for r in full if d87y_positive(r["AMR genotypes core"])]
    positive_target = [r for r in positives if canonical_cluster(r["SNP cluster"]) == TARGET]
    positive_else = [r for r in positives if canonical_cluster(r["SNP cluster"]) != TARGET]
    pct = 100.0 * len(positive_target) / len(positives)
    all_cluster = [r for r in full if canonical_cluster(r["SNP cluster"]) == TARGET]
    print(f"[MASTER] gyrA_D87Y-positive total = {len(positives)}")
    print(f"[MASTER] gyrA_D87Y-positive in {TARGET} = {len(positive_target)}")
    print(f"[MASTER] gyrA_D87Y-positive outside target cluster = {len(positive_else)}")
    print(f"[MASTER] fraction of positives in target cluster = {pct:.6f}%")
    print(f"[MASTER] genomes in target cluster regardless of gyrA_D87Y = {len(all_cluster)}")
    print(f"[COMPARE] claimed 13162/13486: numerator matches={len(positive_target)==13162}, denominator matches={len(positives)==13486}")

    result = {
        "input_master_n": len(full),
        "input_selected_pre_qc_n": len(selected),
        "input_selected_post_qc_n": len(selected_study),
        "target_cluster": TARGET,
        "clade3_n": len(clade3),
        "clade3_target_n": len(target_c3),
        "clade3_target_pct": round(100*len(target_c3)/len(clade3), 6),
        "target_in_other_final_clades": len(target_else),
        "master_gyrA_D87Y_positive_n": len(positives),
        "master_gyrA_D87Y_target_n": len(positive_target),
        "master_gyrA_D87Y_elsewhere_n": len(positive_else),
        "master_fraction_pct": round(pct, 6),
        "master_target_cluster_all_genomes": len(all_cluster),
        "claim_master_numerator_correct": len(positive_target)==13162,
        "claim_master_denominator_correct": len(positives)==13486,
    }
    output = ROOT / "results/clade3_cluster_audit.json"
    output.parent.mkdir(exist_ok=True, parents=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print("[AUDIT]", json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
