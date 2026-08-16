#!/usr/bin/env python3
"""Create compact selected-gene matrices for the Gate10S breast path."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


MATRIX_FILES = [
    "GSE190772_BoM_logCounts.txt.gz",
    "GSM6870693_BoM7_scRNA_LogCounts.txt.gz",
    "GSM6870694_BoM8_scRNA_LogCounts.txt.gz",
]

EXTRA_MARKERS = {
    "LST1", "TYROBP", "FCER1G", "CTSS", "AIF1",
    "CD3D", "CD3E", "TRAC", "CD247",
    "CD4", "IL7R", "LTB", "CD8A", "CD8B",
    "KLRD1", "NCR1", "NCAM1", "FCGR3A",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_markers(project: Path) -> set[str]:
    markers = set(EXTRA_MARKERS)
    signature_path = project / "config/gate9_frozen_ecotype_definition.tsv"
    gate_path = project / "config/gate9a_raw_state_gates.tsv"
    with signature_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            markers.update(row["analysis_markers"].split("|"))
    with gate_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            markers.update(row["core_markers"].split("|"))
    return markers


def extract(source: Path, destination: Path, markers: set[str]) -> dict[str, object]:
    selected: list[tuple[str, list[str]]] = []
    with gzip.open(source, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[-1] != "gene":
            raise ValueError(f"Unexpected header in {source}")
        cells = header[:-1]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(f"Row-width mismatch in {source}")
            gene = fields[-1]
            if gene in markers:
                selected.append((gene, fields[:-1]))
    found = {gene for gene, _ in selected}
    missing = sorted(markers - found)
    with gzip.open(destination, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", *cells])
        for gene, values in selected:
            writer.writerow([gene, *values])
    return {
        "source": source.name,
        "output": destination.name,
        "cells": len(cells),
        "markers_requested": len(markers),
        "markers_found": len(found),
        "duplicate_selected_gene_rows": len(selected) - len(found),
        "missing_markers": "|".join(missing),
        "sha256": sha256(destination),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    project = Path(args.project).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    source_dir = project / "data/raw/GSE190772"
    markers = load_markers(project)
    receipts = []
    for name in MATRIX_FILES:
        source = source_dir / name
        destination = out_dir / name.replace(".txt.gz", ".selected_markers.tsv.gz")
        receipts.append(extract(source, destination, markers))
    receipt_path = out_dir / "selected_marker_extraction_receipt.tsv"
    with receipt_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(receipts[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(receipts)
    outputs = sorted(p for p in out_dir.iterdir() if p.name != "SHA256SUMS")
    with (out_dir / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in outputs:
            handle.write(f"{sha256(path)}  {path.name}\n")
    print("GATE10S_BREAST_MARKER_EXTRACTION=PASS")


if __name__ == "__main__":
    main()
