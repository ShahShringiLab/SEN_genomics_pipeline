#!/usr/bin/env python3
from __future__ import annotations

import csv
import io
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

PHEN = ROOT / "data/part_a/phenotype/FDASEN_AMRPhenotype2.csv"
MASTER = ROOT / "data/part_a/genotype/US-SEN_12192025.csv"
KEYS = ROOT / "data/part_a/selection/Selected_updated_Isolate.csv"
FINAL = ROOT / "data/part_a/selection/Finalupdated_SEN_Genomes1.csv"
QR_CLADE = ROOT / "data/part_a/selection/Selected_US-SEN_12192025_QCAUDIT_R1_qr_clade.csv"
PART_B = ROOT / "metadata/SEN_Genomes.csv"


def require_file(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise AssertionError(f"Missing or empty input: {path}")


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise AssertionError(f"Missing header: {path}")
        rows = list(reader)
        return list(reader.fieldnames), rows


def normalized_master_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    text = path.read_text(encoding="utf-8-sig")
    placeholder = "__SEN_LITERAL_COMMA__"
    if placeholder in text:
        raise AssertionError("Reserved parser placeholder already occurs in master input")

    escaped_comma_n = text.count(r"\,")
    exudate_bad = ',""Exudate, Left Heel"","clinical",'
    malformed_exudate_n = text.count(exudate_bad)

    text = text.replace(r"\,", placeholder)
    text = text.replace(
        exudate_bad,
        ',"Exudate, Left Heel","clinical",',
    )

    reader = csv.reader(io.StringIO(text))
    parsed = list(reader)
    if not parsed:
        raise AssertionError("Master genotype CSV is empty")

    width = len(parsed[0])
    bad = [(i + 1, len(row)) for i, row in enumerate(parsed) if len(row) != width]
    if bad:
        raise AssertionError(f"Master CSV has inconsistent field counts: {bad[:10]}")

    header = parsed[0]
    rows: list[dict[str, str]] = []
    for values in parsed[1:]:
        row = dict(zip(header, values))
        for key, value in row.items():
            row[key] = value.replace(placeholder, ",")
        if "Isolation source" in row and row["Isolation source"].endswith('""'):
            row["Isolation source"] = row["Isolation source"][:-2]
        rows.append(row)

    print(
        "[INFO] master parser normalization:",
        f"{escaped_comma_n} escaped-comma sequence(s);",
        f"{malformed_exudate_n} doubled-quote field(s)",
    )
    return header, rows


def norm(value: str | None) -> str:
    cleaned = " ".join((value or "").split())
    # Match readr handling: literal "NA" and an empty field are both missing.
    return "" if cleaned.upper() == "NA" else cleaned


def unique_nonempty(rows: list[dict[str, str]], field: str) -> set[str]:
    return {row[field] for row in rows if row.get(field, "") != ""}


def main() -> None:
    for path in (PHEN, MASTER, KEYS, FINAL, QR_CLADE, PART_B):
        require_file(path)

    phen_header, phen = read_csv_rows(PHEN)
    master_header, master = normalized_master_rows(MASTER)
    key_header, keys = read_csv_rows(KEYS)
    final_header, final = read_csv_rows(FINAL)
    qr_clade_header, qr_clade = read_csv_rows(QR_CLADE)
    part_b_header, part_b = read_csv_rows(PART_B)

    assert len(phen) == 58, len(phen)
    assert {row["Source"] for row in phen} == {"Human", "Chicken"}
    assert len({(row["Year"], row["Source"]) for row in phen}) == 58
    assert all(row["Cip_Tested"] == row["DSC_Tested"] for row in phen)

    assert len(master_header) == 51, len(master_header)
    assert len(master) == 33831, len(master)
    assert len(unique_nonempty(master, "Isolate")) == 33831
    assert len(unique_nonempty(master, "Run")) == 33831

    assert key_header == ["Strain", "Isolate"], key_header
    assert len(keys) == 3434, len(keys)
    assert len(unique_nonempty(keys, "Isolate")) == 3434
    key_pairs = {(norm(r.get("Strain")), norm(r.get("Isolate"))) for r in keys}
    assert len(key_pairs) == 3434

    assert len(final_header) == 69, len(final_header)
    assert len(final) == 3434, len(final)
    assert len(unique_nonempty(final, "Isolate")) == 3434
    assert len(unique_nonempty(final, "Run")) == 3434
    assert {int(row["QR_genotype"]) for row in final} <= {0, 1}

    assert len(qr_clade) == 3434, len(qr_clade)
    assert len(unique_nonempty(qr_clade, "Run")) == 3434
    assert len(unique_nonempty(qr_clade, "Isolate")) == 3434
    corrected_qr = Counter(int(row["QR_genotype"]) for row in qr_clade)
    assert corrected_qr == Counter({1: 2288, 0: 1146}), corrected_qr
    assert sum(row["QR_corrected"] == "YES" for row in qr_clade) == 3
    audit_clades = Counter(row["Clade"] for row in qr_clade if row["Clade"])
    assert audit_clades == Counter({
        "Clade 1A": 169,
        "Clade 1B": 759,
        "Clade 2": 894,
        "Clade 3": 1357,
        "Unassigned": 127,
    }), audit_clades

    selected = [
        row for row in master
        if (norm(row.get("Strain")), norm(row.get("Isolate"))) in key_pairs
    ]
    assert len(selected) == 3434, len(selected)
    assert len(unique_nonempty(selected, "Isolate")) == 3434

    key_isolates = unique_nonempty(keys, "Isolate")
    final_isolates = unique_nonempty(final, "Isolate")
    assert key_isolates == final_isolates

    assert len(part_b) == 3434, len(part_b)
    assert len(unique_nonempty(part_b, "Run")) == 3434
    assert unique_nonempty(final, "Run") == unique_nonempty(part_b, "Run")

    expected_source = {
        "Human": 1976,
        "Chicken-US": 1059,
        "Other-US": 328,
        "International": 71,
    }
    assert Counter(row["Source"] for row in final) == Counter(expected_source)

    print("[PASS] Part A inputs and the 3,434-genome Part A -> Part B handoff validate.")


if __name__ == "__main__":
    main()
