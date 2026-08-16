#!/usr/bin/env python3
"""Create patient-aware cohort and overlap audits from GEO SOFT metadata."""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "metadata"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def characteristic(text: str, key: str) -> str:
    match = re.search(rf"(?:^|\|\s*){re.escape(key)}:\s*([^|]+)", text, re.I)
    return match.group(1).strip() if match else ""


def parse_patient_samples(samples: list[dict[str, str]]) -> list[dict[str, object]]:
    parsed: list[dict[str, object]] = []
    for row in samples:
        accession = row["accession"]
        title = row["title"]
        patient = ""
        compartment = ""
        disease = ""
        include = "yes"
        exclusion = ""

        if accession == "GSE266330":
            patient = re.sub(r"_[12]$", "", title)
            compartment = "benign" if patient.startswith("ctrl") else "bone_metastasis"
            disease = characteristic(row["characteristics"], "group")
            if not disease:
                disease = patient.split("_", 1)[0]
        elif accession == "GSE143791":
            patient = characteristic(row["characteristics"], "subject id")
            if title.startswith("BMET"):
                compartment = title.split("-", 1)[1].lower()
                disease = "prostate_cancer_bone_metastasis"
            elif title.startswith("BMM"):
                compartment = "benign"
                disease = "benign_bone_marrow"
            else:
                include = "no"
                exclusion = "unrelated add-on human or mouse sample"
                compartment = "other"
                disease = characteristic(row["characteristics"], "phenotype")
        elif accession == "GSE202813":
            patient = characteristic(row["characteristics"], "patient id")
            if title.startswith("RCC-"):
                compartment = title.rsplit("-", 1)[1].lower()
                compartment = {"involve": "involved", "noninvolved": "distal"}.get(
                    compartment, compartment
                )
                disease = "renal_cancer_bone_metastasis"
            else:
                compartment = "benign_stroma" if "stroma" in title else "benign_immune"
                disease = "benign_bone_marrow"
        elif accession == "GSE225209":
            patient = title
            compartment = "bone_metastasis"
            disease = "lung_cancer_bone_metastasis"
        elif accession == "GSE225208":
            patient = title
            if title.startswith("lu-b"):
                compartment = "lung_primary_with_bone_metastasis"
            elif title.startswith("lu-f"):
                compartment = "intrathoracic_lung_cancer"
            else:
                compartment = "bone_metastasis_pool_or_sample"
            disease = "lung_cancer"
        elif accession == "GSE190772":
            patient = title.split()[0] if title.startswith("single cell") else title
            if title.startswith("BoM") or "bone metastasis" in title:
                compartment = "bone_metastasis"
            elif "organoid" in title:
                compartment = "organoid"
            else:
                compartment = "primary"
            disease = "breast_cancer"
        elif accession == "GSE39494":
            patient = title
            compartment = "bone_metastasis" if title.startswith("bone biopsy") else "primary"
            disease = "breast_cancer"
        elif accession == "GSE14776":
            match = re.search(r"Patient\s+(\d+)", title)
            patient = f"Patient{match.group(1)}" if match else title
            compartment = "disseminated_tumor_cell" if title.startswith("DTC") else "metastatic_tumor_cell"
            disease = "breast_cancer"

        parsed.append(
            {
                "accession": accession,
                "gsm": row["gsm"],
                "patient_id_within_study": patient,
                "compartment": compartment,
                "disease": disease,
                "organism": row["organism"],
                "include_primary_analysis": include,
                "exclusion_reason": exclusion,
                "title": title,
            }
        )
    return parsed


def paired_counts(rows: list[dict[str, object]], accession: str) -> tuple[int, int]:
    by_patient: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["accession"] == accession and row["include_primary_analysis"] == "yes":
            by_patient[str(row["patient_id_within_study"])].add(str(row["compartment"]))
    complete = sum({"tumor", "involved", "distal"}.issubset(parts) for parts in by_patient.values())
    any_cancer = sum(any(x in parts for x in {"tumor", "involved", "distal"}) for parts in by_patient.values())
    return complete, any_cancer


def build_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    p143_complete, p143_total = paired_counts(rows, "GSE143791")
    p202_complete, p202_total = paired_counts(rows, "GSE202813")
    return [
        {
            "accession": "GSE143791",
            "role": "paired-gradient discovery/meta-analysis",
            "modality": "scRNA-seq counts",
            "cancer": "prostate",
            "independent_tumor_patients": p143_total,
            "complete_tumor_involved_distal_patients": p143_complete,
            "benign_donors": 7,
            "usable_as_independent_cohort": "yes",
            "key_limitation": "BMM2-BMM9 controls probably overlap GSE202813; two patients have incomplete gradients",
        },
        {
            "accession": "GSE202813",
            "role": "paired-gradient discovery/meta-analysis",
            "modality": "scRNA-seq counts",
            "cancer": "renal",
            "independent_tumor_patients": p202_total,
            "complete_tumor_involved_distal_patients": p202_complete,
            "benign_donors": 9,
            "usable_as_independent_cohort": "yes",
            "key_limitation": "only four complete gradients; benign immune/stroma enrichments are not whole-marrow replicates",
        },
        {
            "accession": "GSE266330",
            "role": "pan-cancer external validation",
            "modality": "scRNA-seq counts plus author objects",
            "cancer": "eight coded cancer groups",
            "independent_tumor_patients": "34 in published design; 42 IDs in updated table",
            "complete_tumor_involved_distal_patients": 0,
            "benign_donors": 5,
            "usable_as_independent_cohort": "conditional",
            "key_limitation": "patient-count conflict and cancer-label inconsistencies require author-ID audit; never infer effects from cell-level replication",
        },
        {
            "accession": "OEP005136",
            "role": "large external pan-cancer replication",
            "modality": "scRNA-seq",
            "cancer": "13 origins",
            "independent_tumor_patients": 52,
            "complete_tumor_involved_distal_patients": "not applicable",
            "benign_donors": "paired/adjacent subset",
            "usable_as_independent_cohort": "conditional on processed-data access",
            "key_limitation": "published in 2026 and already integrates GSE143791/GSE202813; those public samples must be removed before external validation",
        },
        {
            "accession": "GSE225209",
            "role": "lung single-cell validation",
            "modality": "snRNA/scRNA-seq counts",
            "cancer": "lung",
            "independent_tumor_patients": 3,
            "complete_tumor_involved_distal_patients": 0,
            "benign_donors": 0,
            "usable_as_independent_cohort": "validation only",
            "key_limitation": "small n and no normal or primary comparator",
        },
        {
            "accession": "GSE190772",
            "role": "breast single-cell validation",
            "modality": "scRNA-seq plus bulk",
            "cancer": "breast",
            "independent_tumor_patients": "at least 1; BoM7/BoM8 mapping unavailable",
            "complete_tumor_involved_distal_patients": 0,
            "benign_donors": 0,
            "usable_as_independent_cohort": "validation only",
            "key_limitation": "BoM1/BoM2 are bilateral lesions from one patient; GEO does not expose patient mapping for BoM7/BoM8",
        },
        {
            "accession": "GSE225208",
            "role": "lung bulk cross-modal validation",
            "modality": "bulk RNA-seq",
            "cancer": "lung",
            "independent_tumor_patients": "unclear from GEO titles",
            "complete_tumor_involved_distal_patients": 0,
            "benign_donors": 0,
            "usable_as_independent_cohort": "validation only",
            "key_limitation": "28 libraries include likely pooled g-* samples; no pairing should be assumed without paper-level mapping",
        },
        {
            "accession": "GSE39494/GSE14776",
            "role": "legacy breast orthogonal validation",
            "modality": "microarray",
            "cancer": "breast",
            "independent_tumor_patients": "5+11 unique (with within-series overlaps)",
            "complete_tumor_involved_distal_patients": 0,
            "benign_donors": 0,
            "usable_as_independent_cohort": "secondary validation only",
            "key_limitation": "old platform, small sample size, and different isolated-cell/biopsy compositions",
        },
    ]


def main() -> int:
    samples = read_tsv(META / "geo_samples.tsv")
    parsed = parse_patient_samples(samples)
    write_tsv(
        META / "sample_level_audit.tsv",
        parsed,
        [
            "accession",
            "gsm",
            "patient_id_within_study",
            "compartment",
            "disease",
            "organism",
            "include_primary_analysis",
            "exclusion_reason",
            "title",
        ],
    )
    write_tsv(
        META / "cohort_summary.tsv",
        build_summary(parsed),
        [
            "accession",
            "role",
            "modality",
            "cancer",
            "independent_tumor_patients",
            "complete_tumor_involved_distal_patients",
            "benign_donors",
            "usable_as_independent_cohort",
            "key_limitation",
        ],
    )
    overlaps = [
        {
            "source_a": "GSE143791",
            "source_b": "GSE202813",
            "entity": "benign donors BMM2,BMM3,BMM4,BMM5,BMM6,BMM8,BMM9",
            "status": "confirmed shared donors",
            "evidence": "BMM2-BMM9 have thousands of identical normalized 10x cell barcodes across both accessions",
            "required_action": "count each donor once; never use the two accessions as independent normal-control replications; repeat analyses without benign controls",
        },
        {
            "source_a": "OEP005136",
            "source_b": "GSE143791/GSE202813",
            "entity": "49 public paired bone-metastasis/normal samples",
            "status": "documented reuse in 2026 atlas",
            "evidence": "2026 article explicitly lists both GEO accessions as incorporated public data",
            "required_action": "remove reused cells/samples before treating OEP005136 as external validation",
        },
        {
            "source_a": "GSE266330",
            "source_b": "GSE266330 updated clinical table",
            "entity": "patient and cancer labels",
            "status": "internal conflict",
            "evidence": "published design says 34 patients; updated table has 42 tumor IDs and several subtype/origin conflicts",
            "required_action": "use stable patient IDs and author-confirmed labels; run label-agnostic and label-restricted sensitivity analyses",
        },
    ]
    write_tsv(
        META / "overlap_audit.tsv",
        overlaps,
        ["source_a", "source_b", "entity", "status", "evidence", "required_action"],
    )

    counts = Counter((r["accession"], r["include_primary_analysis"]) for r in parsed)
    print(f"COMPLETE samples={len(parsed)} cohort_rows=8 overlap_flags={len(overlaps)}")
    for key, value in sorted(counts.items()):
        print("\t".join((*key, str(value))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
