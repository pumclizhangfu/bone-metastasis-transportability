#!/usr/bin/env python3
"""Audit the official GSE323357 archive without inspecting candidate-axis genes."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Iterable


EXPECTED_SIZE = 974_684_160
SPATIAL = {
    "GSM9564255": "24664",
    "GSM9564256": "24665",
    "GSM9564257": "24666",
    "GSM9564258": "24667",
}
SCRNA = "GSM9564259"


@dataclass
class Check:
    check: str
    status: str
    detail: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_name(name: str) -> bool:
    path = PurePosixPath(name)
    return not path.is_absolute() and ".." not in path.parts


def base_name(name: str) -> str:
    return PurePosixPath(name).name


def open_member_gzip(archive: tarfile.TarFile, member: tarfile.TarInfo):
    raw = archive.extractfile(member)
    if raw is None:
        raise ValueError(f"Cannot read archive member: {member.name}")
    return gzip.GzipFile(fileobj=raw, mode="rb")


def iter_text_lines(handle: BinaryIO) -> Iterable[str]:
    first = True
    for raw_line in handle:
        encoding = "utf-8-sig" if first else "utf-8"
        first = False
        yield raw_line.decode(encoding).rstrip("\r\n")


def gzip_lines(archive: tarfile.TarFile, member: tarfile.TarInfo) -> list[str]:
    with open_member_gzip(archive, member) as handle:
        return list(iter_text_lines(handle))


def matrix_market_dimensions(
    archive: tarfile.TarFile, member: tarfile.TarInfo
) -> tuple[int, int, int]:
    with open_member_gzip(archive, member) as handle:
        lines = iter_text_lines(handle)
        banner = next(lines)
        if not banner.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unexpected Matrix Market banner in {member.name}")
        for line in lines:
            if line and not line.startswith("%"):
                fields = line.split()
                if len(fields) != 3:
                    raise ValueError(f"Malformed Matrix Market dimensions in {member.name}")
                return tuple(map(int, fields))  # type: ignore[return-value]
    raise ValueError(f"Missing Matrix Market dimensions in {member.name}")


def split_row(line: str) -> list[str]:
    if "\t" in line:
        return line.split("\t")
    return next(csv.reader([line]))


def barcode_set(lines: list[str]) -> set[str]:
    return {split_row(line)[0].strip() for line in lines if line.strip()}


def metadata_barcodes(lines: list[str], expected: set[str]) -> tuple[set[str], str]:
    if not lines:
        return set(), "empty"
    rows = [split_row(line) for line in lines if line.strip()]
    header = rows[0]
    body = rows[1:]
    best_index = 0
    best_overlap = -1
    best_values: set[str] = set()
    best_transform = "identity"
    for index in range(max(map(len, rows))):
        raw_values = {row[index].strip() for row in body if len(row) > index}
        candidates = {
            "identity": raw_values,
            "strip_integrated_section_suffix": {
                re.sub(r"_\d+$", "", value) for value in raw_values
            },
        }
        for transform, values in candidates.items():
            overlap = len(values & expected)
            if overlap > best_overlap:
                best_index = index
                best_overlap = overlap
                best_values = values
                best_transform = transform
    label = header[best_index] if len(header) > best_index else f"column_{best_index + 1}"
    return best_values, f"{label};transform={best_transform}"


def position_barcodes(lines: list[str]) -> set[str]:
    values: set[str] = set()
    for line in lines:
        if not line.strip():
            continue
        first = split_row(line)[0].strip()
        if first.lower() not in {"barcode", "barcodes"}:
            values.add(first)
    return values


def consume_gzip(archive: tarfile.TarFile, member: tarfile.TarInfo) -> None:
    with open_member_gzip(archive, member) as handle:
        for _ in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            pass


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tar", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    tar_path = Path(args.tar).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    checks: list[Check] = []
    dimension_rows: list[dict[str, object]] = []
    coordinate_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []

    observed_size = tar_path.stat().st_size
    checks.append(Check("archive_size", "PASS" if observed_size == EXPECTED_SIZE else "FAIL", f"observed={observed_size}; expected={EXPECTED_SIZE}"))
    digest = sha256(tar_path)

    with tarfile.open(tar_path, "r:") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        unsafe = [member.name for member in members if not safe_name(member.name)]
        checks.append(Check("safe_archive_paths", "PASS" if not unsafe else "FAIL", "none" if not unsafe else "|".join(unsafe)))
        by_base = {base_name(member.name): member for member in members}
        accessions = sorted({match.group(0) for name in by_base for match in [re.search(r"GSM\d+", name)] if match})
        expected_accessions = sorted([*SPATIAL, SCRNA])
        checks.append(Check("expected_accessions", "PASS" if accessions == expected_accessions else "FAIL", f"observed={'|'.join(accessions)}; expected={'|'.join(expected_accessions)}"))

        manifest_rows = [
            {"member": member.name, "bytes": member.size, "gzip": str(member.name.endswith(".gz")).upper()}
            for member in members
        ]
        write_tsv(out_dir / "archive_manifest.tsv", manifest_rows, ["member", "bytes", "gzip"])

        for gsm, library in SPATIAL.items():
            prefix = f"{gsm}_Bone_ST_filtered_{library}"
            names = {
                "matrix": f"{prefix}_matrix.mtx.gz",
                "features": f"{prefix}_features.tsv.gz",
                "barcodes": f"{prefix}_barcodes.tsv.gz",
                "positions": f"{gsm}_Bone_ST_{library}_tissue_positions_list.csv.gz",
                "metadata": f"{gsm}_Bone_ST_{library}_spot_metadata.tsv.gz",
                "scalefactors": f"{gsm}_Bone_ST_{library}_scalefactors_json.json.gz",
                "hires": f"{gsm}_Bone_ST_{library}_tissue_hires_image.png.gz",
                "lowres": f"{gsm}_Bone_ST_{library}_tissue_lowres_image.png.gz",
            }
            missing = [name for name in names.values() if name not in by_base]
            checks.append(Check(f"{gsm}_required_files", "PASS" if not missing else "FAIL", "none" if not missing else "|".join(missing)))
            if missing:
                sample_rows.append({"sample": gsm, "modality": "Visium", "library": library, "unit": "section", "status": "FAIL", "detail": "missing required files"})
                continue

            matrix_rows, matrix_cols, nnz = matrix_market_dimensions(archive, by_base[names["matrix"]])
            features = gzip_lines(archive, by_base[names["features"]])
            barcodes = gzip_lines(archive, by_base[names["barcodes"]])
            dimension_status = matrix_rows == len(features) and matrix_cols == len(barcodes)
            dimension_rows.append({"sample": gsm, "modality": "Visium", "matrix_rows": matrix_rows, "feature_rows": len(features), "matrix_cols": matrix_cols, "barcode_rows": len(barcodes), "nnz": nnz, "status": "PASS" if dimension_status else "FAIL"})
            checks.append(Check(f"{gsm}_10x_dimensions", "PASS" if dimension_status else "FAIL", f"matrix={matrix_rows}x{matrix_cols}; features={len(features)}; barcodes={len(barcodes)}; nnz={nnz}"))

            filtered = barcode_set(barcodes)
            positions = position_barcodes(gzip_lines(archive, by_base[names["positions"]]))
            metadata, metadata_column = metadata_barcodes(gzip_lines(archive, by_base[names["metadata"]]), filtered)
            missing_positions = filtered - positions
            missing_metadata = filtered - metadata
            coordinate_status = not missing_positions and not missing_metadata
            coordinate_rows.append({"sample": gsm, "filtered_barcodes": len(filtered), "position_barcodes": len(positions), "metadata_barcodes": len(metadata), "metadata_barcode_column": metadata_column, "missing_in_positions": len(missing_positions), "missing_in_metadata": len(missing_metadata), "status": "PASS" if coordinate_status else "FAIL"})
            checks.append(Check(f"{gsm}_barcode_coverage", "PASS" if coordinate_status else "FAIL", f"filtered={len(filtered)}; positions={len(positions)}; metadata={len(metadata)}; missing_positions={len(missing_positions)}; missing_metadata={len(missing_metadata)}"))

            for key in ("scalefactors", "hires", "lowres"):
                consume_gzip(archive, by_base[names[key]])
            with open_member_gzip(archive, by_base[names["scalefactors"]]) as handle:
                json.load(io.TextIOWrapper(handle, encoding="utf-8"))
            sample_rows.append({"sample": gsm, "modality": "Visium", "library": library, "unit": "section_not_animal", "status": "PASS" if dimension_status and coordinate_status else "FAIL", "detail": "repository animal identifier absent"})

        scrna_names = {
            "matrix": f"{SCRNA}_Bone_scRNA_filtered_counts_matrix.mtx.gz",
            "features": f"{SCRNA}_Bone_scRNA_filtered_counts_features.tsv.gz",
            "barcodes": f"{SCRNA}_Bone_scRNA_filtered_counts_barcodes.tsv.gz",
            "metadata": f"{SCRNA}_Bone_scRNA_cell_metadata.tsv.gz",
        }
        missing = [name for name in scrna_names.values() if name not in by_base]
        checks.append(Check(f"{SCRNA}_required_files", "PASS" if not missing else "FAIL", "none" if not missing else "|".join(missing)))
        if not missing:
            matrix_rows, matrix_cols, nnz = matrix_market_dimensions(archive, by_base[scrna_names["matrix"]])
            features = gzip_lines(archive, by_base[scrna_names["features"]])
            barcodes = gzip_lines(archive, by_base[scrna_names["barcodes"]])
            dimension_status = matrix_rows == len(features) and matrix_cols == len(barcodes)
            dimension_rows.append({"sample": SCRNA, "modality": "scRNA", "matrix_rows": matrix_rows, "feature_rows": len(features), "matrix_cols": matrix_cols, "barcode_rows": len(barcodes), "nnz": nnz, "status": "PASS" if dimension_status else "FAIL"})
            filtered = barcode_set(barcodes)
            metadata, metadata_column = metadata_barcodes(gzip_lines(archive, by_base[scrna_names["metadata"]]), filtered)
            missing_metadata = filtered - metadata
            metadata_status = not missing_metadata
            coordinate_rows.append({"sample": SCRNA, "filtered_barcodes": len(filtered), "position_barcodes": "NA", "metadata_barcodes": len(metadata), "metadata_barcode_column": metadata_column, "missing_in_positions": "NA", "missing_in_metadata": len(missing_metadata), "status": "PASS" if metadata_status else "FAIL"})
            checks.append(Check(f"{SCRNA}_10x_dimensions", "PASS" if dimension_status else "FAIL", f"matrix={matrix_rows}x{matrix_cols}; features={len(features)}; barcodes={len(barcodes)}; nnz={nnz}"))
            checks.append(Check(f"{SCRNA}_barcode_coverage", "PASS" if metadata_status else "FAIL", f"filtered={len(filtered)}; metadata={len(metadata)}; missing_metadata={len(missing_metadata)}"))
            sample_rows.append({"sample": SCRNA, "modality": "scRNA", "library": "scRNA_lib1", "unit": "single_library", "status": "PASS" if dimension_status and metadata_status else "FAIL", "detail": "one repository library"})
        else:
            sample_rows.append({"sample": SCRNA, "modality": "scRNA", "library": "scRNA_lib1", "unit": "single_library", "status": "FAIL", "detail": "missing required files"})

    write_tsv(out_dir / "tenx_dimension_audit.tsv", dimension_rows, ["sample", "modality", "matrix_rows", "feature_rows", "matrix_cols", "barcode_rows", "nnz", "status"])
    write_tsv(out_dir / "barcode_coordinate_audit.tsv", coordinate_rows, ["sample", "filtered_barcodes", "position_barcodes", "metadata_barcodes", "metadata_barcode_column", "missing_in_positions", "missing_in_metadata", "status"])
    write_tsv(out_dir / "sample_unit_audit.tsv", sample_rows, ["sample", "modality", "library", "unit", "status", "detail"])
    check_rows = [{"check": row.check, "status": row.status, "detail": row.detail} for row in checks]
    write_tsv(out_dir / "input_checks.tsv", check_rows, ["check", "status", "detail"])

    overall = "PASS" if checks and all(row.status == "PASS" for row in checks) and all(row["status"] == "PASS" for row in sample_rows) else "FAIL"
    report = [
        "# Gate12E GSE323357 input audit",
        "",
        f"- Overall: **{overall}**",
        f"- Archive: `{tar_path}`",
        f"- Bytes: `{observed_size}` (expected `{EXPECTED_SIZE}`)",
        f"- SHA-256: `{digest}`",
        f"- Archive files: `{len(manifest_rows)}`",
        "- Repository structure: four Visium sections plus one scRNA-seq library.",
        "- Biological-unit boundary: unique animal identifiers are absent; the four sections are descriptive section-level units and cannot be treated as four independent animals.",
        "- Interpretation boundary: mouse spatial support cannot be reported as human spatial validation.",
        "- Candidate-axis expression was not queried during this audit.",
        "",
        "## Checks",
        "",
        "| Check | Status | Detail |",
        "|---|---:|---|",
    ]
    for row in checks:
        report.append(f"| {row.check} | {row.status} | {row.detail.replace('|', ', ')} |")
    report.extend(["", f"`GATE12E_INPUT_AUDIT={overall}`", ""])
    (out_dir / "GATE12E_INPUT_AUDIT.md").write_text("\n".join(report), encoding="utf-8")
    (out_dir / "archive_sha256.txt").write_text(f"{digest}  {tar_path.name}\n", encoding="utf-8")
    print(f"GATE12E_INPUT_AUDIT={overall}")
    print(f"ARCHIVE_SHA256={digest}")


if __name__ == "__main__":
    main()
