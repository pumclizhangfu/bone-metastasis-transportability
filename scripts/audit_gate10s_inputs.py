#!/usr/bin/env python3
"""Audit Gate10S archives without calculating biological scores."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import tarfile
from pathlib import Path


PRIMARY_SIGNATURES = {
    "CD14HI_MONO": "FCN1|VCAN|MNDA|S100A12|CD36|IRAK3|CSF3R|CYBB|S100A8|RBM47|SLC8A1|PLXDC2|CLEC12A|FGD4|AIF1",
    "CD16HI_MONO": "FCGR3A|SMIM25|LST1|FCN1|SERPINA1|TCF7L2|MTSS1|HCK|MS4A7|SLC8A1|CLEC7A|LILRB2|WARS|CLEC12A|IRAK3",
    "MACROPHAGE": "MSR1|DOCK4|SLC8A1|C1QA|CLEC7A|PLXDC2|C1QC|MS4A6A|CD86|HLA-DQA1|FMNL2|MARCKS|C1QB|MAFB|FCGR2A",
    "OSTEOCLAST": "CTSK|ACP5|MMP9|SLC9B2|ATP6V0D2|ANPEP|SPP1|TNFRSF11A|JDP2|CD109|AK5|CSF1R|NRP2|COL27A1|LINC02725",
    "CD4_TREG": "IKZF2|CTLA4|FOXP3|TIGIT|TBC1D4|ICOS|AC013652.1|STAM|BATF|LTB|TRAC|IL32|DUSP16|CD2",
    "CD8_TEX": "TOX|TNFRSF9|DGKH|LINC01934|CD27|CD8A|BICDL1|TNIP3|MIR155HG|NCALD|TTN|GZMK|CD2|HNRNPLL|NKG7",
}

LUNG_SAMPLES = [
    ("GSM7041480", "sg1", "sz", "LUCA_BM_01"),
    ("GSM7041481", "sg2", "s13", "LUCA_BM_02"),
    ("GSM7041482", "sg3", "s14", "LUCA_BM_03"),
]

BREAST_SAMPLES = [
    ("GSM5731346", "BoM1", "BRCA_BM_01", "GSE190772_BoM_logCounts.txt.gz", "GSE190772_BoM_MetaData.txt.gz"),
    ("GSM5731347", "BoM2", "BRCA_BM_01", "GSE190772_BoM_logCounts.txt.gz", "GSE190772_BoM_MetaData.txt.gz"),
    ("GSM6870693", "BoM7", "BRCA_BM_02", "GSM6870693_BoM7_scRNA_LogCounts.txt.gz", "GSM6870693_BoM7_scRNA_MetaData.txt.gz"),
    ("GSM6870694", "BoM8", "BRCA_BM_02", "GSM6870694_BoM8_scRNA_LogCounts.txt.gz", "GSM6870694_BoM8_scRNA_MetaData.txt.gz"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_expected_sha(path: Path) -> str:
    return path.read_text(encoding="utf-8").strip().split()[0]


def gzip_text_from_tar(tar: tarfile.TarFile, member: str):
    extracted = tar.extractfile(member)
    if extracted is None:
        raise ValueError(f"Missing tar member: {member}")
    return io.TextIOWrapper(gzip.GzipFile(fileobj=extracted), encoding="utf-8")


def read_mtx_dims_from_tar(tar: tarfile.TarFile, member: str) -> tuple[int, int, int]:
    with gzip_text_from_tar(tar, member) as handle:
        for line in handle:
            if not line.startswith("%"):
                values = tuple(int(x) for x in line.split())
                if len(values) != 3:
                    raise ValueError(f"Invalid Matrix Market dimensions: {member}")
                return values
    raise ValueError(f"No Matrix Market dimension line: {member}")


def read_features_from_tar(tar: tarfile.TarFile, member: str) -> tuple[int, set[str]]:
    genes: set[str] = set()
    count = 0
    with gzip_text_from_tar(tar, member) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"Invalid features row: {member}")
            count += 1
            genes.add(fields[1])
    return count, genes


def count_barcodes_from_tar(tar: tarfile.TarFile, member: str) -> tuple[int, bool]:
    values = []
    with gzip_text_from_tar(tar, member) as handle:
        for line in handle:
            values.append(line.rstrip("\n"))
    return len(values), len(values) == len(set(values))


def scan_dense_matrix(path: Path) -> tuple[int, int, set[str], list[str], bool]:
    genes: set[str] = set()
    rows = 0
    row_width_ok = True
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if len(header) < 2 or header[-1] != "gene":
            raise ValueError(f"Unexpected dense-matrix header: {path}")
        cells = header[:-1]
        expected = len(header)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            rows += 1
            if len(fields) != expected:
                row_width_ok = False
            if fields:
                genes.add(fields[-1])
    return rows, len(cells), genes, cells, row_width_ok


def scan_metadata(path: Path) -> tuple[int, dict[str, int], set[str]]:
    counts: dict[str, int] = {}
    cell_ids: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing metadata header: {path}")
        origin_col = "orig.ident"
        id_col = "Cell" if "Cell" in reader.fieldnames else "cell.ID"
        for row in reader:
            origin = row[origin_col]
            counts[origin] = counts.get(origin, 0) + 1
            cell_ids.add(row[id_col])
    return sum(counts.values()), counts, cell_ids


def signature_rows(dataset: str, matrix_id: str, genes: set[str]) -> list[dict[str, object]]:
    result = []
    for state, marker_string in PRIMARY_SIGNATURES.items():
        markers = list(dict.fromkeys(marker_string.split("|")))
        present = [gene for gene in markers if gene in genes]
        coverage = len(present) / len(markers)
        result.append({
            "dataset": dataset,
            "matrix_id": matrix_id,
            "state_id": state,
            "markers_expected": len(markers),
            "markers_present": len(present),
            "coverage": coverage,
            "pass_ge_0_80": coverage >= 0.80,
            "missing_markers": "|".join(gene for gene in markers if gene not in genes),
        })
    return result


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    project = Path(args.project).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    lung_dir = project / "data/raw/GSE225209"
    breast_dir = project / "data/raw/GSE190772"
    lung_tar_path = lung_dir / "GSE225209_RAW.tar"
    artifact_rows: list[dict[str, object]] = []
    crosswalk_rows: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []

    expected_tar_sha = read_expected_sha(lung_tar_path.with_suffix(lung_tar_path.suffix + ".sha256"))
    observed_tar_sha = sha256(lung_tar_path)
    artifact_rows.append({
        "dataset": "GSE225209", "artifact": lung_tar_path.name, "check": "sha256",
        "observed": observed_tar_sha, "expected": expected_tar_sha, "pass": observed_tar_sha == expected_tar_sha,
    })

    with tarfile.open(lung_tar_path, "r") as tar:
        members = set(tar.getnames())
        for gsm, stem, sample_name, patient_id in LUNG_SAMPLES:
            prefix = f"{gsm}_{stem}"
            feature_member = f"{prefix}-features.tsv.gz"
            barcode_member = f"{prefix}-barcodes.tsv.gz"
            matrix_member = f"{prefix}-matrix.mtx.gz"
            required = {feature_member, barcode_member, matrix_member}
            complete = required.issubset(members)
            if not complete:
                raise ValueError(f"Incomplete 10X triplet for {sample_name}")
            feature_n, genes = read_features_from_tar(tar, feature_member)
            barcode_n, barcode_unique = count_barcodes_from_tar(tar, barcode_member)
            dims = read_mtx_dims_from_tar(tar, matrix_member)
            dim_pass = dims[0] == feature_n and dims[1] == barcode_n
            artifact_rows.extend([
                {"dataset": "GSE225209", "artifact": prefix, "check": "matrix_dimensions",
                 "observed": "x".join(map(str, dims)), "expected": f"{feature_n}x{barcode_n}xNNZ", "pass": dim_pass},
                {"dataset": "GSE225209", "artifact": prefix, "check": "barcode_uniqueness",
                 "observed": barcode_unique, "expected": True, "pass": barcode_unique},
            ])
            coverage_rows.extend(signature_rows("GSE225209", sample_name, genes))
            crosswalk_rows.append({
                "dataset": "GSE225209", "accession": gsm, "sample_or_lesion": sample_name,
                "matrix_id": prefix, "patient_id": patient_id, "cancer_code": "LUCA",
                "input_type": "raw_10x_counts", "primary_patient_unit": True,
            })

    dense_cache: dict[str, tuple[int, int, set[str], list[str], bool]] = {}
    meta_cache: dict[str, tuple[int, dict[str, int], set[str]]] = {}
    for gsm, lesion, patient_id, matrix_name, metadata_name in BREAST_SAMPLES:
        matrix_path = breast_dir / matrix_name
        metadata_path = breast_dir / metadata_name
        for path in (matrix_path, metadata_path):
            expected = read_expected_sha(path.with_suffix(path.suffix + ".sha256"))
            observed = sha256(path)
            artifact_rows.append({
                "dataset": "GSE190772", "artifact": path.name, "check": "sha256",
                "observed": observed, "expected": expected, "pass": observed == expected,
            })
        if matrix_name not in dense_cache:
            dense_cache[matrix_name] = scan_dense_matrix(matrix_path)
            coverage_rows.extend(signature_rows("GSE190772", matrix_name, dense_cache[matrix_name][2]))
        if metadata_name not in meta_cache:
            meta_cache[metadata_name] = scan_metadata(metadata_path)
        rows, matrix_cells, _, matrix_cell_ids, width_ok = dense_cache[matrix_name]
        meta_rows, origin_counts, metadata_cell_ids = meta_cache[metadata_name]
        exact_cell_set = set(matrix_cell_ids) == metadata_cell_ids
        if matrix_name == "GSE190772_BoM_logCounts.txt.gz":
            lesion_matrix_cells = sum(cell.startswith(f"{lesion}_") for cell in matrix_cell_ids)
            lesion_meta_cells = origin_counts.get(lesion, 0)
        else:
            lesion_matrix_cells = matrix_cells
            lesion_meta_cells = origin_counts.get(lesion, 0)
        artifact_rows.extend([
            {"dataset": "GSE190772", "artifact": f"{lesion}:{matrix_name}", "check": "row_width_consistency",
             "observed": width_ok, "expected": True, "pass": width_ok},
            {"dataset": "GSE190772", "artifact": f"{lesion}:{matrix_name}", "check": "matrix_metadata_cells",
             "observed": lesion_matrix_cells, "expected": lesion_meta_cells, "pass": lesion_matrix_cells == lesion_meta_cells},
            {"dataset": "GSE190772", "artifact": matrix_name, "check": "matrix_metadata_cell_ids",
             "observed": exact_cell_set, "expected": True, "pass": exact_cell_set},
        ])
        crosswalk_rows.append({
            "dataset": "GSE190772", "accession": gsm, "sample_or_lesion": lesion,
            "matrix_id": matrix_name, "patient_id": patient_id, "cancer_code": "BRCA",
            "input_type": "author_log_normalized", "primary_patient_unit": False,
        })

    # De-duplicate repeated whole-file SHA rows caused by two lesions sharing a file.
    artifact_rows = list({(r["dataset"], r["artifact"], r["check"]): r for r in artifact_rows}.values())
    s0_pass = all(bool(r["pass"]) for r in artifact_rows)
    s1_pass = all(bool(r["pass_ge_0_80"]) for r in coverage_rows)
    s2_pass = (
        len([r for r in crosswalk_rows if r["dataset"] == "GSE225209"]) == 3
        and len({r["patient_id"] for r in crosswalk_rows if r["dataset"] == "GSE225209"}) == 3
        and len([r for r in crosswalk_rows if r["dataset"] == "GSE190772"]) == 4
        and len({r["patient_id"] for r in crosswalk_rows if r["dataset"] == "GSE190772"}) == 2
    )
    decision = "READY_FOR_SMOKE" if s0_pass and s1_pass and s2_pass else "STOP_INPUT_AUDIT"

    write_tsv(out_dir / "input_artifact_audit.tsv", artifact_rows)
    write_tsv(out_dir / "sample_lesion_crosswalk.tsv", crosswalk_rows)
    write_tsv(out_dir / "signature_coverage.tsv", coverage_rows)
    receipt = {
        "version": "gate10s_supportive_projection_v1",
        "decision": decision,
        "S0_archive_integrity": s0_pass,
        "S1_signature_coverage": s1_pass,
        "S2_sample_mapping": s2_pass,
        "lung_patients": 3,
        "breast_lesions": 4,
        "breast_patients": 2,
    }
    (out_dir / "input_audit_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    final = [
        "# Gate10S input audit",
        "",
        f"- Decision: **{decision}**",
        f"- S0 archive and dimension integrity: {'PASS' if s0_pass else 'FAIL'}",
        f"- S1 frozen signature coverage >=80%: {'PASS' if s1_pass else 'FAIL'}",
        f"- S2 patient/lesion crosswalk: {'PASS' if s2_pass else 'FAIL'}",
        "- Expression-based biological endpoints calculated: **NO**",
        "- Gate10A decision changed: **NO**",
        "- Gate10B authorized: **NO**",
        "",
        "The audit authorizes only the prespecified smoke tests when READY_FOR_SMOKE.",
        "GSE225209 and GSE190772 remain separate, supportive datasets.",
        "",
    ]
    (out_dir / "GATE10S_INPUT_AUDIT.md").write_text("\n".join(final), encoding="utf-8")
    output_files = sorted(p for p in out_dir.iterdir() if p.name != "SHA256SUMS")
    with (out_dir / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in output_files:
            handle.write(f"{sha256(path)}  {path.name}\n")
    print(f"GATE10S_INPUT_DECISION={decision}")


if __name__ == "__main__":
    main()
