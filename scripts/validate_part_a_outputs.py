#!/usr/bin/env python3
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results/part_a"


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise AssertionError(f"Missing or empty Part A output: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise AssertionError(f"Missing CSV header: {path}")
        return list(reader)


def unique_nonempty(rows: list[dict[str, str]], field: str) -> set[str]:
    return {row.get(field, "") for row in rows if row.get(field, "")}


def as_int(value: str) -> int:
    return int(float(value))


def main() -> None:
    phenotype = read_rows(
        RESULTS / "phenotype/FQ_trends_Human_vs_Chicken.csv"
    )
    if len(phenotype) != 165:
        raise AssertionError(
            f"Phenotype trend output rows={len(phenotype)}; expected 165."
        )
    if {r["Metric"] for r in phenotype} != {
        "Nal-R/DSC or Cip-R", "Nal-R/DSC", "Cip-R"
    }:
        raise AssertionError("Unexpected phenotype metric set.")

    rebuilt = read_rows(
        RESULTS / "cohort/Finalupdated_SEN_Genomes_rebuilt.csv"
    )
    if len(rebuilt) != 3434:
        raise AssertionError(f"Rebuilt cohort rows={len(rebuilt)}; expected 3434.")
    if len(unique_nonempty(rebuilt, "Isolate")) != 3434:
        raise AssertionError("Rebuilt cohort Isolate values are not unique/complete.")
    if len(unique_nonempty(rebuilt, "Run")) != 3434:
        raise AssertionError("Rebuilt cohort Run values are not unique/complete.")

    qr_counts = Counter(as_int(r["QR_genotype"]) for r in rebuilt)
    if qr_counts != Counter({1: 2288, 0: 1146}):
        raise AssertionError(f"Unexpected corrected QR genotype counts: {dict(qr_counts)}")

    missing_clades = {"", "NA"}  # readr writes missing clade assignments as NA
    clade_counts = Counter(
        r.get("Clade", "") for r in rebuilt
        if r.get("Clade", "").strip() not in missing_clades
    )
    excluded_clades = sum(
        r.get("Clade", "").strip() in missing_clades for r in rebuilt
    )
    if excluded_clades != 128:
        raise AssertionError(
            f"Expected 128 excluded genomes; observed {excluded_clades}"
        )
    expected_clades = Counter({
        "Clade 1A": 169,
        "Clade 1B": 759,
        "Clade 2": 894,
        "Clade 3": 1357,
        "Unassigned": 127,
    })
    if clade_counts != expected_clades:
        raise AssertionError(f"Unexpected clade counts: {dict(clade_counts)}")

    summary = read_rows(
        RESULTS / "cohort/cohort_rebuild_validation_summary.csv"
    )
    summary_map = {r["check"]: as_int(r["value"]) for r in summary}
    required_summary = {
        "selected_rows": 3434,
        "unique_isolates": 3434,
        "unique_runs": 3434,
        "collection_year_mismatches": 0,
        "source_mismatches": 0,
        "origin_mismatches": 0,
        "source_group_mismatches": 0,
        "QR_genotype_corrections": 3,
        "QR_marker_positive": 2288,
        "QR_marker_negative": 1146,
        "final_phylogeny_clade_mapped": 3306,
        "final_phylogeny_excluded": 128,
    }
    if summary_map != required_summary:
        raise AssertionError(
            f"Cohort validation summary differs from release contract: {summary_map}"
        )

    marker = read_rows(
        RESULTS / "genotype_summary/FQR_marker_frequency_percentage_95CI_with_susceptible.csv"
    )
    if not marker:
        raise AssertionError("Marker-frequency output is empty.")
    if "Q/FQ-susceptible" not in {r["token_pretty"] for r in marker}:
        raise AssertionError("Marker-frequency output lacks susceptible category.")

    source_rows = read_rows(
        RESULTS / "source_profiles/SEN_Source_counts.csv"
    )
    observed_sources = Counter(
        {r["Source"]: as_int(r["n"]) for r in source_rows}
    )
    expected_sources = Counter({
        "Human": 1976,
        "Chicken-US": 1059,
        "Other-US": 328,
        "International": 71,
    })
    if observed_sources != expected_sources:
        raise AssertionError(
            f"Source counts differ from release contract: {dict(observed_sources)}"
        )

    diversity = read_rows(
        RESULTS / "diversity/QR_diversity_q0_q1_q2_bySource_with95CI.csv"
    )
    expected_div_sources = {"Human", "Chicken-US", "Other-US", "International"}
    if {r["Source"] for r in diversity} != expected_div_sources:
        raise AssertionError("Diversity output does not contain all four source groups.")

    expected_figures = [
        RESULTS / "phenotype/FQ_trends_Human_vs_Chicken.png",
        RESULTS / "phenotype/FQ_trends_Human_vs_Chicken.pdf",
        RESULTS / "genotype_summary/FQR_marker_counts_top100.png",
        RESULTS / "genotype_summary/QRFQR_percent_by_source_year.png",
        RESULTS / "source_profiles/SEN_International_OtherUS_animal_QRFQR.png",
        RESULTS / "diversity/QR_diversity_and_composition.png",
    ]
    for path in expected_figures:
        if not path.is_file() or path.stat().st_size == 0:
            raise AssertionError(f"Missing or empty expected figure: {path}")

    print("[PASS] Part A generated-output contracts validate.")


if __name__ == "__main__":
    main()
