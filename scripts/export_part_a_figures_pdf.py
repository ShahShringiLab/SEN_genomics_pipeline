#!/usr/bin/env python3
"""Export existing figure TIFFs as separate PDF files, without replotting.

PDF output is a lossless image-based rendering of each TIFF (not vector).
An existing native Figure 2 PDF is retained. No original TIFF is changed.
"""
from __future__ import annotations

import sys
from pathlib import Path

try:
    from PIL import Image
except ImportError as exc:
    raise SystemExit(
        "[ERROR] Pillow is needed for TIFF-to-PDF export. Install it in "
        "the sen_stats environment (conda install -n sen_stats pillow), "
        "then rerun this script."
    ) from exc

FIGURES = (
    "FQR_marker_counts_top100",
    "FQ_burden",
    "FQ_burden_5yr",
    "FQ_QRFQR_vs_QFQsusceptible_PERCENT_by_year",
    "FQ_QRFQR_vs_QFQsusceptible_COUNT_by_year",
    "FQ_QSFQS_TOP5_QRFQR_PERCENT_by_year",
    "FQ_QSFQS_TOP5_QRFQR_COUNT_by_year",
)


def export(directory: Path) -> None:
    missing = [name for name in FIGURES if not (directory / (name + ".tif")).is_file()]
    if missing:
        raise SystemExit("[ERROR] Missing original TIFF exports: " + ", ".join(missing))
    for name in FIGURES:
        source = directory / (name + ".tif")
        destination = directory / (name + ".pdf")
        if destination.exists():
            print("[KEEP] Existing PDF: " + str(destination))
            continue
        with Image.open(source) as tif:
            if getattr(tif, "n_frames", 1) != 1:
                raise SystemExit("[ERROR] Multi-frame TIFF unexpected: " + str(source))
            dpi = tif.info.get("dpi", (300, 300))
            try:
                resolution = float(dpi[0])
                if resolution <= 0:
                    resolution = 300.0
            except (ValueError, TypeError, IndexError):
                resolution = 300.0
            if tif.mode == "RGB":
                pdf_image = tif.copy()
            elif tif.mode in ("RGBA", "LA") or "transparency" in tif.info:
                rgba = tif.convert("RGBA")
                pdf_image = Image.new("RGB", rgba.size, "white")
                pdf_image.paste(rgba, mask=rgba.getchannel("A"))
            else:
                pdf_image = tif.convert("RGB")
            pdf_image.save(destination, "PDF", resolution=resolution)
        print("[PDF] " + str(destination))
    if not all((directory / (name + ".pdf")).is_file() for name in FIGURES):
        raise SystemExit("[ERROR] PDF conversion incomplete")
    print("[PASS] Seven separate PDF figures available.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python scripts/export_part_a_figures_pdf.py OUTPUT_DIR")
    export(Path(sys.argv[1]).expanduser().resolve())
