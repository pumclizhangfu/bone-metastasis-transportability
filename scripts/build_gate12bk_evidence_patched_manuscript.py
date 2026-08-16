#!/usr/bin/env python3
"""Apply only Gate12BJ evidence-backed manuscript corrections and audit the delta."""

from __future__ import annotations

import csv
import difflib
import hashlib
import json
import os
import platform
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


COMMAND_PREFIX = "python3 scripts/build_gate12bk_evidence_patched_manuscript.py"

OLD_INTRO = (
    "Spatial transcriptomics adds anatomical context, but spot resolution, deconvolution, cross-species mapping "
    "and uncertain animal replication restrict the claims that can be made <sup>22–26</sup>."
)
DATA_INSERTION_ANCHOR = (
    "Independent single-cell datasets were obtained from GEO under GSE266330 and from the National Omics Data "
    "Encyclopedia project OEP005136."
)
OLD_DATA_PLACEHOLDER = (
    " Exact OEP005136 source-protocol identifiers and the stable publisher download path must be transcribed "
    "from the source record during final journal formatting."
)
OLD_ETHICS_PLACEHOLDER = (
    "Exact protocol identifiers for OEP005136 remain to be transcribed from the final publisher information."
)
GSE323357_IACUC_SENTENCE = (
    "Animal procedures underlying GSE323357 were performed by the source investigators under protocols approved "
    "by the Baylor College of Medicine Institutional Animal Care and Use Committee."
)
AUTHOR_CONTROLLED_TEXT = [
    "**Authors:** [AUTHOR LIST TO BE COMPLETED]",
    "**Affiliations:** [AFFILIATIONS TO BE COMPLETED]",
    "**Corresponding author:** [FULL CONTACT DETAILS TO BE COMPLETED]",
    "[FUNDING SOURCES, GRANT NUMBERS AND FUNDER ROLES TO BE COMPLETED. If no specific funding supported the work, this must be confirmed by all authors.]",
    "[AUTHOR CONFIRMATION REQUIRED. Suggested wording if accurate: “The authors declare no competing interests.”]",
    "[AUTHOR NAMES AND CRediT ROLES TO BE COMPLETED WITHOUT INVENTING CONTRIBUTIONS: Conceptualization, Data curation, Formal analysis, Funding acquisition, Investigation, Methodology, Project administration, Resources, Software, Supervision, Validation, Visualization, Writing – original draft, Writing – review and editing.]",
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def load_config(path: Path) -> dict[str, str]:
    rows = read_tsv(path)
    if not rows or set(rows[0]) != {"key", "value"}:
        raise RuntimeError("Gate12BK configuration must have key and value columns")
    config: dict[str, str] = {}
    for row in rows:
        if row["key"] in config:
            raise RuntimeError(f"Duplicate configuration key: {row['key']}")
        config[row["key"]] = row["value"]
    return config


def expand_citation(token: str) -> list[int]:
    if not re.fullmatch(r"\d+(?:(?:,|–|-)\d+)*", token):
        return []
    numbers: list[int] = []
    for part in token.split(","):
        match = re.fullmatch(r"(\d+)[–-](\d+)", part)
        if match:
            start, end = map(int, match.groups())
            numbers.extend(range(start, end + 1))
        else:
            numbers.append(int(part))
    return numbers


def citation_numbers(text: str) -> list[int]:
    body = text.split("# References\n", 1)[0]
    output: list[int] = []
    for token in re.findall(r"<sup>([^<]+)</sup>", body):
        output.extend(expand_citation(token))
    return output


def reference_records(text: str) -> list[tuple[int, str]]:
    section = text.split("# References\n", 1)[1].split("# Figure legends\n", 1)[0]
    records: list[tuple[int, str]] = []
    for match in re.finditer(r"(?ms)^(\d+)\.\s+(.*?)(?=^\d+\.\s+|\Z)", section.strip()):
        doi_match = re.search(r"https?://doi\.org/([^\s]+)", match.group(2), flags=re.I)
        if not doi_match:
            raise RuntimeError(f"Reference {match.group(1)} has no DOI")
        records.append((int(match.group(1)), doi_match.group(1).rstrip(".").lower()))
    return records


def section_between(text: str, start: str | None, end: str | None) -> str:
    start_index = 0 if start is None else text.index(start)
    end_index = len(text) if end is None else text.index(end, start_index + (len(start) if start else 0))
    return text[start_index:end_index]


def add_check(
    rows: list[dict[str, Any]],
    check_id: str,
    domain: str,
    observed: Any,
    expected: Any,
    passed: bool,
    note: str = "",
) -> None:
    rows.append(
        {
            "check_id": check_id,
            "domain": domain,
            "observed": observed,
            "expected": expected,
            "status": "PASS" if passed else "FAIL",
            "note": note,
        }
    )


def output_manifest(root: Path, excluded: set[str]) -> list[dict[str, Any]]:
    rows = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative in excluded:
            continue
        rows.append({"relative_path": relative, "bytes": path.stat().st_size, "sha256": sha256_file(path)})
    return rows


def main() -> int:
    started = datetime.now().astimezone()
    if len(sys.argv) != 2:
        raise RuntimeError(f"Usage: {COMMAND_PREFIX} <config.tsv>")
    root = Path.cwd().resolve()
    config_path = Path(sys.argv[1]).resolve()
    try:
        config_display = config_path.relative_to(root).as_posix()
    except ValueError:
        config_display = str(config_path)
    run_command = f"{COMMAND_PREFIX} {config_display}"
    config = load_config(config_path)
    required_keys = {
        "source_manuscript", "source_legends", "source_figure_manifest", "gate12bj_receipt",
        "gate12bj_correction_plan", "gate12bj_report", "gate12bj_provenance_audit", "gate12bi_blockers",
        "audit_script", "output_dir", "expected_manuscript_sha256", "expected_legends_sha256",
        "expected_figure_manifest_sha256", "expected_gate12bj_receipt_sha256",
        "expected_gate12bj_correction_plan_sha256", "expected_gate12bj_report_sha256",
        "expected_gate12bj_provenance_audit_sha256", "expected_gate12bi_blockers_sha256",
        "expected_audit_script_sha256", "expected_references", "expected_citation_occurrences",
        "expected_figure_legends", "expected_resolved_blockers", "expected_remaining_blockers",
    }
    missing = sorted(required_keys - set(config))
    if missing:
        raise RuntimeError("Missing Gate12BK configuration keys: " + ", ".join(missing))

    paths = {
        "source_manuscript": root / config["source_manuscript"],
        "source_legends": root / config["source_legends"],
        "source_figure_manifest": root / config["source_figure_manifest"],
        "gate12bj_receipt": root / config["gate12bj_receipt"],
        "gate12bj_correction_plan": root / config["gate12bj_correction_plan"],
        "gate12bj_report": root / config["gate12bj_report"],
        "gate12bj_provenance_audit": root / config["gate12bj_provenance_audit"],
        "gate12bi_blockers": root / config["gate12bi_blockers"],
        "audit_script": root / config["audit_script"],
    }
    expected_hashes = {
        "source_manuscript": config["expected_manuscript_sha256"],
        "source_legends": config["expected_legends_sha256"],
        "source_figure_manifest": config["expected_figure_manifest_sha256"],
        "gate12bj_receipt": config["expected_gate12bj_receipt_sha256"],
        "gate12bj_correction_plan": config["expected_gate12bj_correction_plan_sha256"],
        "gate12bj_report": config["expected_gate12bj_report_sha256"],
        "gate12bj_provenance_audit": config["expected_gate12bj_provenance_audit_sha256"],
        "gate12bi_blockers": config["expected_gate12bi_blockers_sha256"],
        "audit_script": config["expected_audit_script_sha256"],
    }
    for role, expected in expected_hashes.items():
        path = paths[role]
        if not path.is_file():
            raise RuntimeError(f"Missing frozen input: {path}")
        observed = sha256_file(path)
        if observed != expected:
            raise RuntimeError(f"Frozen input hash mismatch for {path}: expected {expected}; observed {observed}")

    output = root / config["output_dir"]
    staging = output.parent / f".{output.name}.staging"
    if output.exists():
        raise RuntimeError(f"Refusing to overwrite existing Gate12BK output: {output}")
    if staging.exists():
        raise RuntimeError(f"Refusing to overwrite existing Gate12BK staging directory: {staging}")

    bj_receipt = json.loads(paths["gate12bj_receipt"].read_text(encoding="utf-8"))
    if bj_receipt.get("status") != "COMPLETED":
        raise RuntimeError("Gate12BJ parent receipt is not COMPLETED")
    if bj_receipt.get("verdicts", {}).get("overall") != "PASS_WITH_REQUIRED_CORRECTIONS":
        raise RuntimeError("Gate12BJ parent verdict is not PASS_WITH_REQUIRED_CORRECTIONS")
    if bj_receipt.get("verdicts", {}).get("reference_metadata") != "PASS":
        raise RuntimeError("Gate12BJ reference metadata did not pass")
    if bj_receipt.get("verdicts", {}).get("dataset_provenance") != "PASS":
        raise RuntimeError("Gate12BJ dataset provenance did not pass")

    correction_rows = read_tsv(paths["gate12bj_correction_plan"])
    corrections = {row["correction_id"]: row for row in correction_rows}
    if set(corrections) != {"BJ-01", "BJ-02", "BJ-03", "BJ-04"}:
        raise RuntimeError(f"Unexpected Gate12BJ correction set: {sorted(corrections)}")
    for correction_id in ["BJ-01", "BJ-02", "BJ-03"]:
        if corrections[correction_id]["status"] != "READY_TO_PATCH":
            raise RuntimeError(f"{correction_id} is not READY_TO_PATCH")
    if corrections["BJ-04"]["status"] != "NO_TEXT_CHANGE_REQUIRED":
        raise RuntimeError("BJ-04 is not marked NO_TEXT_CHANGE_REQUIRED")

    source = paths["source_manuscript"].read_text(encoding="utf-8")
    legends = paths["source_legends"].read_text(encoding="utf-8")
    new_intro = corrections["BJ-03"]["recommended_text"]
    new_data = corrections["BJ-01"]["recommended_text"]
    new_ethics = corrections["BJ-02"]["recommended_text"]
    exact_old_fragments = [OLD_INTRO, DATA_INSERTION_ANCHOR, OLD_DATA_PLACEHOLDER, OLD_ETHICS_PLACEHOLDER]
    for fragment in exact_old_fragments:
        count = source.count(fragment)
        if count != 1:
            raise RuntimeError(f"Expected exactly one source occurrence; observed {count}: {fragment[:90]}")
    if source.count(GSE323357_IACUC_SENTENCE) != 1:
        raise RuntimeError("Frozen manuscript does not contain the exact GSE323357 IACUC sentence")

    revised = source.replace(OLD_INTRO, new_intro, 1)
    revised = revised.replace(DATA_INSERTION_ANCHOR, DATA_INSERTION_ANCHOR + " " + new_data, 1)
    revised = revised.replace(OLD_DATA_PLACEHOLDER, "", 1)
    revised = revised.replace(OLD_ETHICS_PLACEHOLDER, new_ethics, 1)

    source_masked = source.replace(OLD_INTRO, "[[BJ-03]]", 1)
    source_masked = source_masked.replace(DATA_INSERTION_ANCHOR, DATA_INSERTION_ANCHOR + " [[BJ-01]]", 1)
    source_masked = source_masked.replace(OLD_DATA_PLACEHOLDER, "", 1)
    source_masked = source_masked.replace(OLD_ETHICS_PLACEHOLDER, "[[BJ-02]]", 1)
    revised_masked = revised.replace(new_intro, "[[BJ-03]]", 1)
    revised_masked = revised_masked.replace(new_data, "[[BJ-01]]", 1)
    revised_masked = revised_masked.replace(new_ethics, "[[BJ-02]]", 1)
    allowed_delta_only = source_masked == revised_masked
    if not allowed_delta_only:
        raise RuntimeError("Revised manuscript contains changes outside BJ-01 to BJ-03")

    staging.mkdir(parents=True)
    manuscript_dir = staging / "manuscript"
    provenance_dir = staging / "provenance"
    manuscript_dir.mkdir()
    provenance_dir.mkdir()
    revised_path = manuscript_dir / "Gate12BK_Evidence_Patched_Manuscript.md"
    legends_path = manuscript_dir / "Gate12BK_Figure_Legends.md"
    revised_path.write_text(revised, encoding="utf-8")
    legends_path.write_bytes(paths["source_legends"].read_bytes())

    change_audit = [
        {
            "correction_id": "BJ-01", "target": "Data availability", "action": "APPLIED",
            "old_occurrences": source.count(OLD_DATA_PLACEHOLDER), "new_occurrences": revised.count(new_data),
            "source_evidence": corrections["BJ-01"]["evidence"],
            "status": "PASS" if source.count(OLD_DATA_PLACEHOLDER) == 1 and revised.count(new_data) == 1 and OLD_DATA_PLACEHOLDER not in revised else "FAIL",
        },
        {
            "correction_id": "BJ-02", "target": "Ethics statement", "action": "APPLIED",
            "old_occurrences": source.count(OLD_ETHICS_PLACEHOLDER), "new_occurrences": revised.count(new_ethics),
            "source_evidence": corrections["BJ-02"]["evidence"],
            "status": "PASS" if source.count(OLD_ETHICS_PLACEHOLDER) == 1 and revised.count(new_ethics) == 1 and OLD_ETHICS_PLACEHOLDER not in revised else "FAIL",
        },
        {
            "correction_id": "BJ-03", "target": "Introduction spatial-limitations sentence", "action": "APPLIED",
            "old_occurrences": source.count(OLD_INTRO), "new_occurrences": revised.count(new_intro),
            "source_evidence": corrections["BJ-03"]["evidence"],
            "status": "PASS" if source.count(OLD_INTRO) == 1 and revised.count(new_intro) == 1 and OLD_INTRO not in revised else "FAIL",
        },
        {
            "correction_id": "BJ-04", "target": "GSE323357 animal-ethics wording", "action": "RETAINED_UNCHANGED",
            "old_occurrences": source.count(GSE323357_IACUC_SENTENCE),
            "new_occurrences": revised.count(GSE323357_IACUC_SENTENCE),
            "source_evidence": corrections["BJ-04"]["evidence"],
            "status": "PASS" if source.count(GSE323357_IACUC_SENTENCE) == revised.count(GSE323357_IACUC_SENTENCE) == 1 else "FAIL",
        },
    ]
    write_tsv(
        staging / "GATE12BK_CHANGE_AUDIT.tsv", change_audit,
        ["correction_id", "target", "action", "old_occurrences", "new_occurrences", "source_evidence", "status"],
    )

    section_specs = [
        ("Front matter", None, "# Abstract\n", "UNCHANGED"),
        ("Abstract", "# Abstract\n", "# Introduction\n", "UNCHANGED"),
        ("Introduction", "# Introduction\n", "# Results\n", "EXPECTED_CHANGE"),
        ("Results", "# Results\n", "# Discussion\n", "UNCHANGED"),
        ("Discussion", "# Discussion\n", "# Conclusions\n", "UNCHANGED"),
        ("Conclusions", "# Conclusions\n", "# Materials and methods\n", "UNCHANGED"),
        ("Materials and methods", "# Materials and methods\n", "# Data availability\n", "UNCHANGED"),
        ("Data availability", "# Data availability\n", "# Ethics statement\n", "EXPECTED_CHANGE"),
        ("Ethics statement", "# Ethics statement\n", "# Acknowledgements and funding\n", "EXPECTED_CHANGE"),
        ("Author-controlled declarations", "# Acknowledgements and funding\n", "# References\n", "UNCHANGED"),
        ("References", "# References\n", "# Figure legends\n", "UNCHANGED"),
        ("Figure legends", "# Figure legends\n", None, "UNCHANGED"),
    ]
    section_rows = []
    for section_name, start, end, expectation in section_specs:
        before = section_between(source, start, end)
        after = section_between(revised, start, end)
        same = before == after
        passed = same if expectation == "UNCHANGED" else not same
        section_rows.append(
            {
                "section": section_name,
                "expected": expectation,
                "source_sha256": sha256_text(before),
                "revised_sha256": sha256_text(after),
                "observed": "UNCHANGED" if same else "CHANGED",
                "status": "PASS" if passed else "FAIL",
            }
        )
    write_tsv(
        staging / "GATE12BK_SECTION_HASH_AUDIT.tsv", section_rows,
        ["section", "expected", "source_sha256", "revised_sha256", "observed", "status"],
    )

    source_citations = citation_numbers(source)
    revised_citations = citation_numbers(revised)
    source_references = reference_records(source)
    revised_references = reference_records(revised)
    expected_references = int(config["expected_references"])
    expected_citations = int(config["expected_citation_occurrences"])
    checks: list[dict[str, Any]] = []
    add_check(checks, "BK-01", "scope", allowed_delta_only, True, allowed_delta_only, "Masked manuscripts must be byte-identical")
    add_check(checks, "BK-02", "corrections", sum(row["action"] == "APPLIED" and row["status"] == "PASS" for row in change_audit), 3, sum(row["action"] == "APPLIED" and row["status"] == "PASS" for row in change_audit) == 3)
    add_check(checks, "BK-03", "corrections", change_audit[3]["status"], "PASS", change_audit[3]["status"] == "PASS", "BJ-04 must remain unchanged")
    add_check(checks, "BK-04", "sections", sum(row["status"] == "PASS" for row in section_rows), len(section_rows), all(row["status"] == "PASS" for row in section_rows))
    add_check(checks, "BK-05", "references", len(revised_references), expected_references, len(revised_references) == expected_references)
    add_check(checks, "BK-06", "references", revised_references, source_references, revised_references == source_references, "Reference order and DOI list must be byte-content invariant")
    add_check(checks, "BK-07", "citations", len(revised_citations), expected_citations, len(revised_citations) == expected_citations)
    add_check(checks, "BK-08", "citations", revised_citations, source_citations, revised_citations == source_citations)
    add_check(checks, "BK-09", "citations", sorted(set(revised_citations)), list(range(1, expected_references + 1)), sorted(set(revised_citations)) == list(range(1, expected_references + 1)))
    add_check(checks, "BK-10", "legends", sha256_file(legends_path), config["expected_legends_sha256"], sha256_file(legends_path) == config["expected_legends_sha256"])
    legend_count = len(re.findall(r"(?m)^## (?:Figure [1-6]|Supplementary Figure S(?:[1-9]|1[0-2]))\.", legends))
    add_check(checks, "BK-11", "legends", legend_count, int(config["expected_figure_legends"]), legend_count == int(config["expected_figure_legends"]))
    add_check(checks, "BK-12", "provenance", revised.count("https://www.biosino.org/node/project/detail/OEP005136"), 1, revised.count("https://www.biosino.org/node/project/detail/OEP005136") == 1)
    add_check(checks, "BK-13", "provenance", revised.count("OEZ00021715"), 1, revised.count("OEZ00021715") == 1)
    add_check(checks, "BK-14", "provenance", revised.count("OED01122886"), 1, revised.count("OED01122886") == 1)
    add_check(checks, "BK-15", "ethics", revised.count("050432-4-2108∗"), 1, revised.count("050432-4-2108∗") == 1)
    add_check(checks, "BK-16", "placeholders", OLD_DATA_PLACEHOLDER in revised or OLD_ETHICS_PLACEHOLDER in revised, False, OLD_DATA_PLACEHOLDER not in revised and OLD_ETHICS_PLACEHOLDER not in revised)
    author_text_preserved = all(source.count(value) == revised.count(value) == 1 for value in AUTHOR_CONTROLLED_TEXT)
    add_check(checks, "BK-17", "author-controlled", author_text_preserved, True, author_text_preserved, "Author metadata and declarations must remain unresolved, not inferred")
    write_tsv(
        staging / "GATE12BK_LOCAL_CHECKS.tsv", checks,
        ["check_id", "domain", "observed", "expected", "status", "note"],
    )
    if any(row["status"] != "PASS" for row in change_audit + section_rows + checks):
        failures = [row for row in change_audit + section_rows + checks if row.get("status") != "PASS"]
        raise RuntimeError(f"Gate12BK local audit failed: {failures}")

    resolved_ids = set(config["expected_resolved_blockers"].split("|"))
    blocker_rows = read_tsv(paths["gate12bi_blockers"])
    disposition_rows = []
    resolution_evidence = {
        "BI-08": "Gate12BK Data availability now supplies the stable OEP005136 URL plus OEZ00021715/OED01122886; Gate12BJ provenance PASS.",
        "BI-09": "Gate12BK Ethics statement now supplies Fudan University Shanghai Cancer Center IRB 050432-4-2108∗ and consent wording; Gate12BJ full-text provenance PASS.",
        "BI-16": "Gate12BJ run_v2 verified 42/42 DOI metadata records and all citation locations.",
    }
    for row in blocker_rows:
        blocker_id = row["blocker_id"]
        disposition_rows.append(
            {
                **row,
                "gate12bk_status": "RESOLVED" if blocker_id in resolved_ids else "UNRESOLVED",
                "gate12bk_evidence": resolution_evidence.get(blocker_id, "Outside Gate12BK evidence-patch scope."),
            }
        )
    unresolved = [row for row in disposition_rows if row["gate12bk_status"] == "UNRESOLVED"]
    if len(unresolved) != int(config["expected_remaining_blockers"]):
        raise RuntimeError(f"Expected {config['expected_remaining_blockers']} remaining blockers; observed {len(unresolved)}")
    write_tsv(
        staging / "GATE12BK_BLOCKER_DISPOSITION.tsv", disposition_rows,
        ["blocker_id", "category", "severity", "finding", "evidence", "required_action", "resolved", "gate12bk_status", "gate12bk_evidence"],
    )

    diff_lines = difflib.unified_diff(
        source.splitlines(keepends=True), revised.splitlines(keepends=True),
        fromfile=config["source_manuscript"],
        tofile=(Path(config["output_dir"]) / "manuscript/Gate12BK_Evidence_Patched_Manuscript.md").as_posix(),
    )
    (staging / "GATE12BK_UNIFIED_DIFF.patch").write_text("".join(diff_lines), encoding="utf-8")

    input_rows = []
    for role, path in paths.items():
        input_rows.append(
            {"input_role": role, "relative_path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
    input_rows.append(
        {"input_role": "config", "relative_path": config_path.relative_to(root).as_posix(), "bytes": config_path.stat().st_size, "sha256": sha256_file(config_path)}
    )
    write_tsv(provenance_dir / "GATE12BK_INPUT_HASHES.tsv", input_rows, ["input_role", "relative_path", "bytes", "sha256"])
    (provenance_dir / "GATE12BK_ENVIRONMENT.txt").write_text(
        "\n".join(
            [
                f"python={sys.version.replace(os.linesep, ' ')}",
                f"platform={platform.platform()}",
                f"working_directory={root}",
                f"command={run_command}",
            ]
        ) + "\n",
        encoding="utf-8",
    )

    unresolved_blocking = sum(row["severity"] == "BLOCKING" for row in unresolved)
    unresolved_required = sum(row["severity"] == "REQUIRED_FIX" for row in unresolved)
    unresolved_target = sum(row["severity"] == "TARGET_DEPENDENT" for row in unresolved)
    report = [
        "# Gate12BK evidence-patched manuscript and local integrity audit",
        "",
        "## Executive verdict",
        "",
        "- Evidence-backed patch integrity: **PASS**.",
        "- Gate12BJ corrections: **BJ-01 to BJ-03 APPLIED; BJ-04 VERIFIED UNCHANGED**.",
        "- Reference/citation integrity: **PASS** (42 references; 47 citation occurrences; no order or DOI change).",
        "- Scientific results, discussion, conclusions, methods and figure legends: **BYTE-CONTENT UNCHANGED**.",
        "- Submission readiness: **NOT_READY**.",
        f"- Remaining Gate12BI items: {len(unresolved)} ({unresolved_blocking} blocking, {unresolved_required} required fixes, {unresolved_target} target-dependent).",
        "",
        "## Applied corrections",
        "",
        "1. Added the stable OEP005136 project URL and the paired NODE object identifiers OEZ00021715/OED01122886.",
        "2. Added the source-reported Fudan University Shanghai Cancer Center IRB 050432-4-2108∗ and informed-consent wording.",
        "3. Separated general spatial-method limitations supported by references 22-26 from this study's cross-species and unavailable-animal-identifier constraints.",
        "4. Retained the verified Baylor IACUC statement without fabricating an unreported protocol number.",
        "",
        "## Resolved Gate12BI items",
        "",
        "- BI-08: OEP005136 source path and object identifiers.",
        "- BI-09: OEP005136 ethics and consent wording.",
        "- BI-16: external reference and citation-support verification.",
        "",
        "## Deliberately unresolved",
        "",
        "Author list, affiliations, corresponding-author details, funding, conflicts and CRediT roles remain untouched because they require author input. Public code/source-data deposition, stale package wording, release-content claims, formatted Supplementary Tables S4/S6 and target-journal formatting/policy checks also remain unresolved. Gate12BK must not be represented as a submission-ready package.",
        "",
        "## Next gate",
        "",
        "Gate12BL should address the non-author-controlled package defects BI-10 to BI-13: correct the stale release wording, reconcile claims about scripts/configurations/intermediate objects with the actual release, and create journal-ready Supplementary Tables S4 and S6. Author-controlled items remain deferred until supplied.",
    ]
    (staging / "GATE12BK_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    ended = datetime.now().astimezone()
    (staging / "GATE12BK_RUN_LOG.txt").write_text(
        "\n".join(
            [
                f"started_at={started.isoformat()}", f"command={run_command}",
                f"source_manuscript_sha256={expected_hashes['source_manuscript']}",
                f"gate12bj_receipt_sha256={expected_hashes['gate12bj_receipt']}",
                "patch_integrity=PASS", "submission_readiness=NOT_READY",
                f"ended_at={ended.isoformat()}",
            ]
        ) + "\n",
        encoding="utf-8",
    )

    manifest_path = provenance_dir / "GATE12BK_OUTPUT_MANIFEST.tsv"
    manifest = output_manifest(staging, {"GATE12BK_RECEIPT.json", "provenance/GATE12BK_OUTPUT_MANIFEST.tsv"})
    write_tsv(manifest_path, manifest, ["relative_path", "bytes", "sha256"])
    receipt = {
        "schema_version": "1.0",
        "gate_id": "Gate12BK",
        "run_id": "gate12bk_run_v1",
        "status": "COMPLETED",
        "command": run_command,
        "working_directory": str(root),
        "started_at": started.isoformat(),
        "ended_at": ended.isoformat(),
        "exit_code": 0,
        "parents": {
            "gate12bh_manuscript_sha256": expected_hashes["source_manuscript"],
            "gate12bj_receipt_sha256": expected_hashes["gate12bj_receipt"],
            "gate12bi_blockers_sha256": expected_hashes["gate12bi_blockers"],
        },
        "verdicts": {
            "evidence_patch_integrity": "PASS",
            "citation_reference_integrity": "PASS",
            "scientific_core_invariance": "PASS",
            "submission_readiness": "NOT_READY",
        },
        "metrics": {
            "applied_corrections": 3,
            "verified_unchanged_corrections": 1,
            "references": len(revised_references),
            "citation_occurrences": len(revised_citations),
            "local_checks_passed": sum(row["status"] == "PASS" for row in checks),
            "local_checks_total": len(checks),
            "resolved_gate12bi_items": sorted(resolved_ids),
            "remaining_gate12bi_items": len(unresolved),
        },
        "output_manuscript": {
            "path": "manuscript/Gate12BK_Evidence_Patched_Manuscript.md",
            "sha256": sha256_file(revised_path),
        },
        "output_manifest": {
            "path": "provenance/GATE12BK_OUTPUT_MANIFEST.tsv",
            "sha256": sha256_file(manifest_path),
            "rows": len(manifest),
        },
    }
    (staging / "GATE12BK_RECEIPT.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    staging.rename(output)
    print(json.dumps({"gate": "Gate12BK", "patch_integrity": "PASS", "submission_readiness": "NOT_READY", "output": str(output)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
