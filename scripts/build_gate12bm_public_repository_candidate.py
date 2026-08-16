#!/usr/bin/env python3
"""Build the Gate12BM sanitized public-repository candidate.

This builder performs no network operation and does not create a remote
repository, choose a license, insert author metadata, or mint a DOI.
"""

from __future__ import annotations

import argparse
import ast
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


SCHEMA_VERSION = "gate12bm-1.0"
EXPECTED_VERDICT = "CANDIDATE_BUILT_WITH_BLOCKERS"
MAX_REPOSITORY_FILE_BYTES = 5 * 1024 * 1024


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def read_key_value_tsv(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"key", "value"}:
        raise ValueError(f"Expected key/value TSV: {path}")
    values: dict[str, str] = {}
    for row in rows:
        key = row["key"].strip()
        if not key or key in values:
            raise ValueError(f"Blank or duplicate config key in {path}: {key!r}")
        values[key] = row["value"].strip()
    return values


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def write_text(path: Path, text: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")
    path.chmod(0o755 if executable else 0o644)


def ensure_relative(path_text: str, label: str) -> Path:
    path = Path(path_text)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"{label} must be a project-relative path: {path_text}")
    return path


def assert_hash(path: Path, expected: str, label: str) -> None:
    observed = sha256_file(path)
    if observed != expected:
        raise RuntimeError(f"{label} hash mismatch: expected {expected}, observed {observed}")


def validate_allowlist(rows: list[dict[str, str]], project: Path, kind: str) -> None:
    required = {"source_path", "destination_path", "stage", "public_status"}
    if not rows or set(rows[0]) != required:
        raise ValueError(f"Invalid {kind} allow-list columns")
    sources: set[str] = set()
    destinations: set[str] = set()
    for row in rows:
        source = ensure_relative(row["source_path"], f"{kind} source")
        destination = ensure_relative(row["destination_path"], f"{kind} destination")
        if row["source_path"] in sources or row["destination_path"] in destinations:
            raise ValueError(f"Duplicate {kind} allow-list entry: {row}")
        sources.add(row["source_path"])
        destinations.add(row["destination_path"])
        if not (project / source).is_file():
            raise FileNotFoundError(f"Missing allow-listed {kind}: {source}")
        if not str(destination).startswith(("scripts/", "config/")):
            raise ValueError(f"Unexpected {kind} destination: {destination}")


def allowlist_content_digest(
    rows: list[dict[str, str]],
    project: Path,
    excluded_sources: set[str] | None = None,
) -> str:
    """Hash the ordered source-path/file-hash pairs behind an allow-list."""
    excluded = excluded_sources or set()
    digest = hashlib.sha256()
    for row in sorted(rows, key=lambda item: item["source_path"]):
        source_path = row["source_path"]
        if source_path in excluded:
            continue
        line = f"{source_path}\t{sha256_file(project / source_path)}\n"
        digest.update(line.encode("utf-8"))
    return digest.hexdigest()


def verify_gate12bl_manifest(bl: Path, manifest: Path) -> int:
    """Verify the complete upstream manifest and exact selected-payload closure."""
    rows = read_tsv(manifest)
    if not rows or set(rows[0]) != {"relative_path", "bytes", "sha256"}:
        raise ValueError("Invalid Gate12BL output manifest")
    listed: set[str] = set()
    for row in rows:
        relative = ensure_relative(row["relative_path"], "Gate12BL manifest path").as_posix()
        if relative in listed:
            raise RuntimeError(f"Duplicate Gate12BL manifest path: {relative}")
        listed.add(relative)
        path = bl / relative
        if not path.is_file():
            raise RuntimeError(f"Gate12BL manifest file missing: {relative}")
        if path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"Gate12BL manifest size mismatch: {relative}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Gate12BL manifest hash mismatch: {relative}")

    selected_prefixes = ("figures/", "source_data/", "tables/", "manuscript/")
    selected_listed = {path for path in listed if path.startswith(selected_prefixes)}
    selected_actual = {
        path.relative_to(bl).as_posix()
        for directory in ("figures", "source_data", "tables", "manuscript")
        for path in (bl / directory).rglob("*")
        if path.is_file()
    }
    if selected_listed != selected_actual:
        missing = sorted(selected_actual - selected_listed)
        stale = sorted(selected_listed - selected_actual)
        raise RuntimeError(f"Gate12BL selected payload is not manifest-closed; missing={missing}, stale={stale}")
    return len(rows)


def sanitize_public_text(text: str) -> tuple[str, list[str]]:
    """Remove private defaults while preserving parseable public code."""
    changes: list[str] = []

    project_pattern = re.compile(r'"/(?:Users|home)/[^"\n]*/multi_cohort_gate2"')
    text, count = project_pattern.subn('"."', text)
    if count:
        changes.append(f"project_root_default_to_current_directory:{count}")

    python_pattern = re.compile(r'"/(?:Users|home)/[^"\n]*/python3"')
    text, count = python_pattern.subn('"python3"', text)
    if count:
        changes.append(f"private_python_runtime_to_python3:{count}")

    office_pattern = re.compile(r'"/(?:Users|home)/[^"\n]*/soffice"')
    text, count = office_pattern.subn('"soffice"', text)
    if count:
        changes.append(f"private_office_runtime_to_soffice:{count}")

    var_prefix = r'"/' + r'var/tmp/[^"\n]*'
    replacements = [
        (var_prefix + r'gate12b"', '"results/gate12b_cell_states"', "gate12b_temp_default"),
        (var_prefix + r'nichenet_prior"', '"data/external/nichenet_prior"', "nichenet_temp_default"),
        (var_prefix + r'raw_matrices"', '"data/raw/gse266330/raw_matrices"', "gse266330_temp_default"),
        (var_prefix + r'mbone_extracted"', '"data/raw/oep005136/mbone_extracted"', "oep_temp_default"),
    ]
    for pattern, replacement, label in replacements:
        text, count = re.subn(pattern, replacement, text)
        if count:
            changes.append(f"{label}:{count}")

    return text, changes


def copy_sanitized_allowlist(
    project: Path,
    repository: Path,
    rows: list[dict[str, str]],
    ledger: list[dict[str, object]],
) -> None:
    for row in rows:
        source = project / row["source_path"]
        destination = repository / row["destination_path"]
        raw = source.read_text(encoding="utf-8")
        public, changes = sanitize_public_text(raw)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if changes:
            destination.write_text(public, encoding="utf-8")
        else:
            shutil.copyfile(source, destination)
        destination.chmod(0o755 if source.suffix == ".sh" else 0o644)
        source_hash = sha256_file(source)
        candidate_hash = sha256_file(destination)
        byte_identical = source_hash == candidate_hash
        if row["public_status"] == "included_exact" and not byte_identical:
            raise RuntimeError(f"Exact allow-list entry was altered: {source}")
        ledger.append(
            {
                "source_path": row["source_path"],
                "candidate_path": f"repository/{row['destination_path']}",
                "stage": row["stage"],
                "declared_status": row["public_status"],
                "transformations": ";".join(changes) if changes else "none",
                "source_sha256": source_hash,
                "candidate_sha256": candidate_hash,
                "byte_identical": str(byte_identical).upper(),
            }
        )


def copy_tree_exact(source: Path, destination: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if not source.is_dir():
        raise FileNotFoundError(source)
    for path in sorted(p for p in source.rglob("*") if p.is_file()):
        relative = path.relative_to(source)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        observed = sha256_file(target)
        expected = sha256_file(path)
        if observed != expected:
            raise RuntimeError(f"Exact-copy failure: {path}")
        rows.append(
            {
                "source_path": str(path),
                "candidate_path": str(target),
                "bytes": target.stat().st_size,
                "sha256": observed,
            }
        )
    return rows


def manifest_rows(root: Path, exclude: set[Path] | None = None) -> list[dict[str, object]]:
    excluded = {p.resolve() for p in (exclude or set())}
    rows: list[dict[str, object]] = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        if path.resolve() in excluded:
            continue
        rows.append(
            {
                "relative_path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return rows


def verify_manifest(root: Path, manifest: Path) -> None:
    rows = read_tsv(manifest)
    for row in rows:
        path = root / row["relative_path"]
        if not path.is_file():
            raise RuntimeError(f"Manifest file missing: {path}")
        if path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"Manifest size mismatch: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Manifest hash mismatch: {path}")


def make_deterministic_tar_gz(source: Path, archive: Path, root_name: str) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    with archive.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as gzip_handle:
            with tarfile.open(fileobj=gzip_handle, mode="w") as tar_handle:
                for path in sorted(p for p in source.rglob("*") if p.is_file()):
                    relative = path.relative_to(source).as_posix()
                    info = tarfile.TarInfo(name=f"{root_name}/{relative}")
                    info.size = path.stat().st_size
                    info.mtime = 0
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    relative_path = path.relative_to(source).as_posix()
                    info.mode = 0o755 if relative_path in {"workflow/run.sh", "tests/accept_release.py"} else 0o644
                    with path.open("rb") as handle:
                        tar_handle.addfile(info, handle)


def audit_tar_members(archive: Path, expected_root: str) -> None:
    seen: set[str] = set()
    with tarfile.open(archive, mode="r:gz") as handle:
        for member in handle.getmembers():
            member_path = Path(member.name)
            if member.name in seen:
                raise RuntimeError(f"Duplicate archive member: {member.name}")
            seen.add(member.name)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise RuntimeError(f"Unsafe archive member: {member.name}")
            if not member_path.parts or member_path.parts[0] != expected_root:
                raise RuntimeError(f"Unexpected archive root: {member.name}")
            if not member.isfile():
                raise RuntimeError(f"Non-regular archive member: {member.name}")


def archive_roundtrip_verify(github_archive: Path, companion_archive: Path) -> str:
    """Extract both distributable archives and run the bundled verifier."""
    audit_tar_members(github_archive, "repository")
    audit_tar_members(companion_archive, "companion")
    with tempfile.TemporaryDirectory(prefix="gate12bm-archive-roundtrip-") as tmp:
        root = Path(tmp)
        with tarfile.open(github_archive, mode="r:gz") as handle:
            for member in handle.getmembers():
                target = root / member.name
                target.parent.mkdir(parents=True, exist_ok=True)
                source = handle.extractfile(member)
                if source is None:
                    raise RuntimeError(f"Could not read archive member: {member.name}")
                with source, target.open("wb") as destination:
                    shutil.copyfileobj(source, destination)
                target.chmod(member.mode)
        with tarfile.open(companion_archive, mode="r:gz") as handle:
            for member in handle.getmembers():
                target = root / member.name
                target.parent.mkdir(parents=True, exist_ok=True)
                source = handle.extractfile(member)
                if source is None:
                    raise RuntimeError(f"Could not read archive member: {member.name}")
                with source, target.open("wb") as destination:
                    shutil.copyfileobj(source, destination)
                target.chmod(member.mode)
        command = [str(root / "repository" / "workflow" / "run.sh"), "verify"]
        result = subprocess.run(command, cwd=root / "repository", text=True, capture_output=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"Archive round-trip verification failed: {result.stderr or result.stdout}")
        return result.stdout + result.stderr


def syntax_audit(repository: Path, rscript: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted((repository / "scripts").rglob("*")):
        if not path.is_file():
            continue
        status = "PASS"
        detail = ""
        try:
            if path.suffix == ".py":
                ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
                detail = "python_ast_parse"
            elif path.suffix.lower() == ".r":
                command = [
                    rscript,
                    "--vanilla",
                    "-e",
                    "parse(file=commandArgs(trailingOnly=TRUE)[1])",
                    str(path),
                ]
                result = subprocess.run(command, text=True, capture_output=True, check=False)
                if result.returncode != 0:
                    raise RuntimeError(result.stderr.strip() or result.stdout.strip())
                detail = "R_parse"
            elif path.suffix == ".sh":
                result = subprocess.run(["bash", "-n", str(path)], text=True, capture_output=True, check=False)
                if result.returncode != 0:
                    raise RuntimeError(result.stderr.strip() or result.stdout.strip())
                detail = "bash_n"
            else:
                continue
        except Exception as exc:  # noqa: BLE001 - audit must record the exact failure
            status = "FAIL"
            detail = str(exc).replace("\t", " ").replace("\n", " ")[:1000]
        rows.append(
            {
                "relative_path": path.relative_to(repository).as_posix(),
                "language": path.suffix.lower().lstrip("."),
                "status": status,
                "detail": detail,
            }
        )
    if any(row["status"] != "PASS" for row in rows):
        raise RuntimeError("One or more public scripts failed syntax parsing")
    return rows


def security_audit(repository: Path, companion: Path) -> list[dict[str, object]]:
    findings: list[dict[str, object]] = []
    forbidden_extensions = {
        ".rds", ".rdata", ".rda", ".h5ad", ".h5", ".bam", ".cram",
        ".fastq", ".fq", ".pyc", ".soft", ".xlsx", ".pdf", ".png",
    }
    literal_rules = {
        "private_macos_path": "/" + "Users/",
        "private_linux_home": "/" + "home/",
        "private_var_tmp": "/" + "var/tmp/",
        "private_storage": "/" + "OceanStor100D/",
        "private_cloud_host": "xiyou" + "cloud",
        "private_webshell": "web" + "shell.",
        "private_sftp_host": "fms." + "biosino.org",
        "ssh_private_key": "BEGIN OPENSSH" + " PRIVATE KEY",
        "pem_private_key": "BEGIN" + " PRIVATE KEY",
        "codex_attachment": ".codex/" + "attachments",
    }
    personal_name = "lizhang" + "fu"
    email_re = re.compile(r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.(?:com|org|edu|gov|net|cn|io)\b", re.IGNORECASE)
    secret_re = re.compile(
        r"(?i)\b(password|passwd|api[_-]?key|access[_-]?token|secret)\b\s*[:=]\s*['\"][^'\"]{8,}"
        r"|\b(?:ghp|github_pat|sk|xox[baprs])[-_A-Za-z0-9]{16,}\b"
    )
    private_ip_re = re.compile(r"\b(?:10(?:\.\d{1,3}){3}|192\.168(?:\.\d{1,3}){2}|172\.(?:1[6-9]|2\d|3[01])(?:\.\d{1,3}){2})\b")

    for path in sorted(repository.rglob("*")):
        if path.is_symlink():
            findings.append({"scope": "repository", "path": str(path.relative_to(repository)), "rule": "symlink", "status": "FAIL", "detail": "symlink not permitted"})
            continue
        if not path.is_file():
            continue
        relative = path.relative_to(repository).as_posix()
        if path.suffix.lower() in forbidden_extensions:
            findings.append({"scope": "repository", "path": relative, "rule": "file_type", "status": "FAIL", "detail": path.suffix.lower()})
        if path.stat().st_size > MAX_REPOSITORY_FILE_BYTES:
            findings.append({"scope": "repository", "path": relative, "rule": "file_size", "status": "FAIL", "detail": str(path.stat().st_size)})
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            findings.append({"scope": "repository", "path": relative, "rule": "non_text_repository_file", "status": "FAIL", "detail": "UTF-8 decode failed"})
            continue
        lowered = text.lower()
        if personal_name in lowered:
            findings.append({"scope": "repository", "path": relative, "rule": "personal_identifier", "status": "FAIL", "detail": "personal user/name token"})
        for label, literal in literal_rules.items():
            if literal.lower() in lowered:
                findings.append({"scope": "repository", "path": relative, "rule": label, "status": "FAIL", "detail": "private literal detected"})
        for email in email_re.findall(text):
            if not email.lower().endswith("@example.org"):
                findings.append({"scope": "repository", "path": relative, "rule": "email", "status": "FAIL", "detail": "non-placeholder email detected"})
        if secret_re.search(text):
            findings.append({"scope": "repository", "path": relative, "rule": "credential_assignment", "status": "FAIL", "detail": "credential-like assignment detected"})
        if private_ip_re.search(text):
            findings.append({"scope": "repository", "path": relative, "rule": "private_ip", "status": "FAIL", "detail": "RFC1918 address detected"})

    disallowed_companion = {".rds", ".rdata", ".rda", ".h5ad", ".h5", ".bam", ".cram", ".fastq", ".fq", ".pyc", ".soft"}
    sensitive_columns = {"age", "sex", "gender", "race", "ethnicity", "date_of_birth", "birth_date", "diagnosis_date", "treatment_date", "death_date", "address", "phone", "email", "treatment", "response", "vital_status", "survival", "overall_survival"}
    for path in sorted(companion.rglob("*")):
        if path.is_symlink():
            findings.append({"scope": "companion", "path": str(path.relative_to(companion)), "rule": "symlink", "status": "FAIL", "detail": "symlink not permitted"})
            continue
        if not path.is_file():
            continue
        relative = path.relative_to(companion).as_posix()
        suffixes = {suffix.lower() for suffix in path.suffixes}
        if suffixes & disallowed_companion:
            findings.append({"scope": "companion", "path": relative, "rule": "restricted_object_type", "status": "FAIL", "detail": ",".join(sorted(suffixes & disallowed_companion))})
        if "/source_data/" in f"/{relative}" and (path.suffix == ".tsv" or path.name.endswith(".tsv.gz")):
            opener = gzip.open if path.name.endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8") as handle:
                header = next(csv.reader(handle, delimiter="\t"), [])
            hits = sorted({column.strip().lower() for column in header} & sensitive_columns)
            if hits:
                findings.append({"scope": "companion", "path": relative, "rule": "sensitive_clinical_column", "status": "FAIL", "detail": ",".join(hits)})

        text_payloads: list[bytes] = []
        try:
            if path.name.endswith(".tsv.gz"):
                with gzip.open(path, "rb") as handle:
                    text_payloads.append(handle.read())
            elif path.suffix.lower() in {".tsv", ".md", ".txt", ".json", ".pdf"}:
                text_payloads.append(path.read_bytes())
            elif path.suffix.lower() == ".xlsx":
                with zipfile.ZipFile(path) as archive:
                    for member in archive.namelist():
                        if member.lower().endswith((".xml", ".rels")):
                            text_payloads.append(archive.read(member))
        except (OSError, EOFError, zipfile.BadZipFile) as exc:
            findings.append({"scope": "companion", "path": relative, "rule": "embedded_text_scan", "status": "FAIL", "detail": str(exc)})
        for payload in text_payloads:
            decoded = payload.decode("latin-1", errors="ignore")
            lowered = decoded.lower()
            if personal_name in lowered:
                findings.append({"scope": "companion", "path": relative, "rule": "personal_identifier", "status": "FAIL", "detail": "personal user/name token"})
            for label, literal in literal_rules.items():
                if literal.lower() in lowered:
                    findings.append({"scope": "companion", "path": relative, "rule": label, "status": "FAIL", "detail": "private literal detected"})
            for email in email_re.findall(decoded):
                if not email.lower().endswith("@example.org"):
                    findings.append({"scope": "companion", "path": relative, "rule": "email", "status": "FAIL", "detail": "non-placeholder email detected"})
            if secret_re.search(decoded):
                findings.append({"scope": "companion", "path": relative, "rule": "credential_assignment", "status": "FAIL", "detail": "credential-like assignment detected"})
            if private_ip_re.search(decoded):
                findings.append({"scope": "companion", "path": relative, "rule": "private_ip", "status": "FAIL", "detail": "RFC1918 address detected"})

    if findings:
        sample = json.dumps(findings[:10], ensure_ascii=False, sort_keys=True)
        raise RuntimeError(f"Security/privacy audit failed; first findings: {sample}")
    pass_rules = [
        ("repository", "symlink", "No symlink"),
        ("repository", "file_type", "No prohibited binary/raw file type"),
        ("repository", "file_size", "No file exceeds 5 MiB"),
        ("repository", "utf8_text", "All repository payloads decode as UTF-8 text"),
        ("repository", "personal_identifier", "No personal user/name token"),
        ("repository", "private_absolute_path", "No private macOS/Linux/tmp/storage path"),
        ("repository", "private_infrastructure", "No private cloud host token"),
        ("repository", "email", "No non-placeholder email"),
        ("repository", "credential_assignment", "No credential-like assignment"),
        ("repository", "private_ip", "No RFC1918 address"),
        ("companion", "symlink", "No symlink"),
        ("companion", "restricted_object_type", "No restricted raw/processed object type"),
        ("companion", "sensitive_clinical_column", "No prohibited clinical/quasi-identifying column"),
        ("companion", "embedded_text_scan", "TSV/GZIP/XLSX XML/Markdown/JSON text and raw PDF bytes scanned"),
        ("companion", "personal_identifier", "No personal user/name token"),
        ("companion", "private_absolute_path", "No private macOS/Linux/tmp/storage path"),
        ("companion", "private_infrastructure", "No private cloud host token"),
        ("companion", "email", "No non-placeholder email"),
        ("companion", "credential_assignment", "No credential-like assignment"),
        ("companion", "private_ip", "No RFC1918 address"),
    ]
    return [
        {"scope": scope, "path": ".", "rule": rule, "status": "PASS", "detail": detail}
        for scope, rule, detail in pass_rules
    ]


def portable_verifier_source() -> str:
    return r'''#!/usr/bin/env python3
"""Offline verifier for the Gate12BM public candidate; standard library only."""
from __future__ import annotations
import csv
import hashlib
import sys
from pathlib import Path

def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def verify(root: Path, manifest_name: str) -> tuple[int, list[str]]:
    manifest = root / "manifests" / manifest_name
    failures = []
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    listed = set()
    for row in rows:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or relative.as_posix() in listed:
            failures.append(f"unsafe_or_duplicate_manifest_path:{row['relative_path']}")
            continue
        listed.add(relative.as_posix())
        path = root / relative
        if not path.is_file():
            failures.append(f"missing:{path}")
        elif path.stat().st_size != int(row["bytes"]):
            failures.append(f"size:{path}")
        elif digest(path) != row["sha256"]:
            failures.append(f"hash:{path}")
    for path in root.rglob("*"):
        if path.is_symlink():
            failures.append(f"symlink:{path}")
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != manifest
    }
    for relative in sorted(actual - listed):
        failures.append(f"unregistered:{root / relative}")
    for relative in sorted(listed - actual):
        failures.append(f"manifest_only:{root / relative}")
    return len(rows), failures

def main() -> int:
    repository = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    companion = repository.parent / "companion"
    repo_n, repo_fail = verify(repository, "repository_files.tsv")
    companion_n, companion_fail = verify(companion, "companion_files.tsv")
    failures = repo_fail + companion_fail
    figure_n = len([p for p in (companion / "reference" / "figures").rglob("*") if p.is_file()])
    source_n = len([p for p in (companion / "data" / "source_data").rglob("*") if p.is_file()])
    if figure_n != 42:
        failures.append(f"figure_count:{figure_n}")
    if source_n != 85:
        failures.append(f"source_data_count:{source_n}")
    print(f"repository_manifest_files={repo_n}")
    print(f"companion_manifest_files={companion_n}")
    print(f"figure_files={figure_n}")
    print(f"source_data_files={source_n}")
    print("verdict=" + ("PASS" if not failures else "FAIL"))
    for item in failures:
        print(item, file=sys.stderr)
    return 0 if not failures else 1

if __name__ == "__main__":
    raise SystemExit(main())
'''


def workflow_source() -> str:
    return r'''#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)"
MODE="${1:-verify}"
case "$MODE" in
  verify)
    exec python3 "$ROOT_DIR/tests/accept_release.py" "$ROOT_DIR"
    ;;
  render|analysis-frozen|full-raw)
    printf '%s\n' "$MODE: DOCUMENTED_NOT_YET_VALIDATED in Gate12BM" >&2
    exit 3
    ;;
  *)
    printf '%s\n' "usage: workflow/run.sh {verify|render|analysis-frozen|full-raw}" >&2
    exit 2
    ;;
esac
'''


def build_repository_documents(repository: Path) -> None:
    write_text(
        repository / "README.md",
        """# Pan-cancer bone-metastasis single-cell and spatial transcriptomics

This directory is a **local deposition candidate**, not an already published
repository. It contains sanitized analysis/rendering code, frozen configuration,
provenance ledgers, and an offline verifier. The adjacent `companion/` directory
contains the accepted static figures and source-data payload.

## Verified in Gate12BM

- frozen Gate12BL figures and source-data files are present and hash-verified;
- allow-listed public scripts parse successfully;
- manifests can be verified offline with Python's standard library;
- private paths, personal identifiers, credentials, disallowed raw objects, and
  repository files larger than 5 MiB are rejected.

Run from this directory:

```bash
workflow/run.sh verify
```

When using the two candidate archives, extract both into the same parent
directory; they create sibling `repository/` and `companion/` directories.

## Not yet validated

`render`, `analysis-frozen`, and `full-raw` intentionally return
`DOCUMENTED_NOT_YET_VALIDATED`. The project currently lacks a validated public
lockfile/container, distributable discovery SCE objects, distributable NicheNet
priors, and an end-to-end clean-environment reconstruction receipt.

## Author actions before public deposit

1. resolve the manuscript disclosure gap for GSE225209 and GSE190772;
2. choose a software/content license after checking third-party constraints;
3. approve redistribution of the companion payload;
4. replace author-controlled citation placeholders;
5. create the remote repository and archive release, then insert the real URL/DOI.

No URL, DOI, license, author, email, ORCID, or affiliation is fabricated here.
""",
    )
    write_text(
        repository / "CHANGELOG.md",
        """# Changelog

## Gate12BM candidate v1

- split code-oriented repository and static-evidence companion payload;
- added frozen script/config allow-lists and source-to-candidate hashes;
- added accession, authoritative-run, method-to-code, and artifact-DAG ledgers;
- added portable offline manifest verification;
- recorded BI-18 for incomplete disclosure of supportive communication datasets;
- left license, redistribution approval, remote URL, and DOI under author control.
""",
    )
    write_text(
        repository / "CITATION.cff.in",
        """cff-version: 1.2.0
message: "If you use this software, cite the associated article and archived release."
title: "Pan-cancer bone-metastasis single-cell and spatial transcriptomics reproducibility materials"
type: software
authors:
  - family-names: "AUTHOR_TO_COMPLETE"
    given-names: "AUTHOR_TO_COMPLETE"
version: "AUTHOR_TO_COMPLETE"
date-released: "AUTHOR_TO_COMPLETE"
repository-code: "AUTHOR_TO_COMPLETE"
doi: "AUTHOR_TO_COMPLETE"
""",
    )
    write_text(
        repository / ".gitignore",
        """.DS_Store
__pycache__/
*.pyc
*.RData
*.Rhistory
*.rds
*.h5ad
data/raw/*
!data/raw/README.md
""",
    )
    write_text(
        repository / "LICENSES" / "LICENSE_DECISION_REQUIRED.md",
        """# License decision required

Gate12BM does not grant a license. Before public deposit, the authors must choose
an appropriate license for original code/content and confirm that third-party
inputs, derived reference objects, figures, and source data may be redistributed
under the intended terms. Do not rename this file to `LICENSE` until that review
is complete.
""",
    )
    write_text(
        repository / "LICENSES" / "THIRD_PARTY_NOTICES.md",
        """# Third-party notices and boundaries

Raw GEO/NODE datasets, source publications, software packages, NicheNet priors,
immune reference objects, and histology/source imagery remain governed by their
respective providers and licenses. Gate12BM does not redistribute restricted
expression matrices, raw sequencing archives, third-party RDS objects, or source
histology assets. Public accession identifiers and transformation provenance are
provided instead.
""",
    )
    write_text(
        repository / "docs" / "REPRODUCIBILITY_SCOPE.md",
        """# Reproducibility scope

The validated level is **static publication-payload verification**: accepted
figures, source-data tables, supplementary tables, and manuscript files are
present and independently rehashed. Script syntax is also checked.

Figure rerendering, processed-input reconstruction, clean-environment execution,
and raw-to-final reconstruction are separate higher levels and remain
`NOT_VALIDATED`. A passing `verify` command must not be cited as evidence that
those higher levels have passed.
""",
    )
    write_text(
        repository / "docs" / "DATA_ACCESS.md",
        """# Data access

Public GEO records can be obtained from the accession URLs listed in
`data/accessions.tsv`. OEP005136/NODE material follows provider-controlled access.
The paired mBone archive used in the study is identified by analysis OEZ00021715
and data object OED01122886; credentials and private SFTP paths are intentionally
excluded.

GSE225209 is a supportive communication co-detection dataset. GSE190772 was
prespecified for the same support layer but was non-evaluable under the primary
eligibility rule. Both must be disclosed in the final manuscript even though the
latter did not contribute an evaluable patient-level estimate.
""",
    )
    write_text(
        repository / "docs" / "KNOWN_LIMITATIONS.md",
        """# Known limitations

- no validated `renv.lock`, Python lockfile, or container image is available;
- development receipts span Ubuntu/R 4.5.1, macOS/R 4.6.1, and Python 3.12.13;
- discovery annotated SCE objects and the final Harmony checkpoint are not in the
  redistributable candidate;
- NicheNet priors and immune-reference RDS objects are third-party inputs and are
  not redistributed;
- NODE/OEP raw matrices require provider-authorized access;
- exact cross-platform graphics hashes are not promised;
- the manuscript currently omits GSE225209/GSE190772 from Data availability;
- no license, remote repository URL, archive DOI, or author metadata is assigned.
""",
    )
    write_text(
        repository / "docs" / "DEPOSITION_CHECKLIST.md",
        """# Deposition checklist

- [ ] Repair BI-18 in the manuscript and re-audit dataset-role disclosure.
- [ ] Confirm author list, affiliations, corresponding author, funding, COI, and CRediT roles.
- [ ] Choose and approve the repository/content license.
- [ ] Review companion files for redistribution authority.
- [ ] Create a clean-environment lock/container and validate at least `verify` plus intended rerender targets.
- [ ] Create the author-controlled GitHub repository.
- [ ] Deposit the frozen companion archive and obtain the real DOI.
- [ ] Replace `CITATION.cff.in` placeholders and the manuscript repository placeholder.
- [ ] Select the journal and complete target-specific formatting/compliance checks.
""",
    )
    write_text(
        repository / "environment" / "ENVIRONMENT_HISTORY.md",
        """# Environment history

Historical computation used Ubuntu 22.04 with R 4.5.1 for server-side analysis,
macOS with R 4.6.1 for later rendering/audits, and Python 3.12.13 for the Gate12BL
release audit. These observations are provenance, not a portable environment
lock. Package versions and system libraries must be frozen and clean-tested in a
future environment gate before raw-to-final reproducibility can be claimed.
""",
    )
    write_text(
        repository / "environment" / "LOCKFILES_NOT_YET_VALIDATED.md",
        """# Lockfiles not yet validated

No `renv.lock`, Python lockfile, Conda environment, or container recipe is
asserted in Gate12BM. Creating an untested lockfile from the current workstation
would give false confidence; a future gate must resolve packages in a clean
environment, run the intended workflow, and retain its receipt.
""",
    )
    write_text(
        repository / "environment" / "system-requirements.txt",
        """Documented historical requirements; not a validated lock:
- R/Bioconductor and Seurat/Harmony/data.table/ggplot2 family packages
- Python 3.12 standard library for candidate verification
- system graphics libraries and fonts for publication rendering
- provider-authorized access for restricted NODE/OEP inputs
""",
    )
    write_text(repository / "environment" / "fonts" / "README.md", "No font binaries are redistributed. Publication rendering used locally available sans-serif fonts; future containers must declare and test their font stack.")
    write_text(repository / "data" / "raw" / "README.md", "Raw data are not redistributed. See `../accessions.tsv` and `../../docs/DATA_ACCESS.md`.")
    write_text(repository / "data" / "frozen_inputs" / "README.md", "Large or restricted processed inputs are not included. Their absence is why `analysis-frozen` and `full-raw` are not validated.")
    write_text(repository / "tests" / "accept_release.py", portable_verifier_source(), executable=True)
    write_text(repository / "workflow" / "run.sh", workflow_source(), executable=True)


def build_ledgers(repository: Path) -> None:
    accession_rows = [
        {"accession_or_id": "GSE143791", "provider": "GEO", "role": "discovery single-cell cohort", "access": "public", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143791", "redistributed": "no"},
        {"accession_or_id": "GSE202813", "provider": "GEO", "role": "external single-cell cohort", "access": "public", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202813", "redistributed": "no"},
        {"accession_or_id": "GSE266330", "provider": "GEO", "role": "external paired/clinical single-cell support", "access": "public", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266330", "redistributed": "no"},
        {"accession_or_id": "GSE323357", "provider": "GEO", "role": "spatial transcriptomic cohort", "access": "public", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE323357", "redistributed": "no"},
        {"accession_or_id": "OEP005136", "provider": "NODE/BMDC", "role": "external bone-metastasis cohort", "access": "provider controlled", "url": "https://www.biosino.org/node/project/detail/OEP005136", "redistributed": "no"},
        {"accession_or_id": "OEZ00021715;OED01122886", "provider": "NODE/BMDC", "role": "paired mBone analysis and archive identifiers", "access": "provider controlled", "url": "https://www.biosino.org/node/project/detail/OEP005136", "redistributed": "no"},
        {"accession_or_id": "GSE225209", "provider": "GEO", "role": "supportive communication co-detection", "access": "public", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225209", "redistributed": "no"},
        {"accession_or_id": "GSE190772", "provider": "GEO", "role": "prespecified communication support; non-evaluable under primary rule", "access": "public", "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190772", "redistributed": "no"},
    ]
    write_tsv(repository / "data" / "accessions.tsv", ["accession_or_id", "provider", "role", "access", "url", "redistributed"], accession_rows)

    run_rows = [
        {"stage": "publication_payload", "authoritative_run": "Gate12BL run_v1", "path_role": "companion static payload", "validation_level": "accepted_hash_verified"},
        {"stage": "main_figures_1", "authoritative_run": "Gate12AZ", "path_role": "Figure 1 renderer", "validation_level": "historical_authoritative"},
        {"stage": "main_figures_2_3", "authoritative_run": "Gate12AZ", "path_role": "Figures 2-3 renderer", "validation_level": "historical_authoritative"},
        {"stage": "main_figures_4_6", "authoritative_run": "Gate12BE", "path_role": "Figures 4-6 renderer", "validation_level": "historical_authoritative"},
        {"stage": "supplementary_base", "authoritative_run": "Gate12AF/Gate12BD", "path_role": "supplementary renderers", "validation_level": "historical_authoritative"},
        {"stage": "supplementary_s2", "authoritative_run": "Gate12BF v3 -> Gate12BG v4", "path_role": "analysis-unit and CLR repair", "validation_level": "historical_authoritative"},
        {"stage": "integrated_submission", "authoritative_run": "Gate12BH", "path_role": "accepted figure/source-data assembly", "validation_level": "accepted"},
        {"stage": "active_manuscript", "authoritative_run": "Gate12BK -> Gate12BL", "path_role": "evidence-patched active manuscript", "validation_level": "accepted_with_open_blockers"},
    ]
    write_tsv(repository / "config" / "authoritative_runs.tsv", ["stage", "authoritative_run", "path_role", "validation_level"], run_rows)

    stage_rows = [
        {"order": 1, "stage": "acquisition", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=acquisition", "status": "documented_not_validated"},
        {"order": 2, "stage": "discovery_atlas", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=discovery", "status": "documented_not_validated"},
        {"order": 3, "stage": "lineage_umap", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=umap", "status": "documented_not_validated"},
        {"order": 4, "stage": "communication", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=communication", "status": "documented_not_validated"},
        {"order": 5, "stage": "supportive_projection", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=supportive_projection", "status": "documented_not_validated"},
        {"order": 6, "stage": "external_validation", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=external_validation", "status": "documented_not_validated"},
        {"order": 7, "stage": "spatial", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=spatial", "status": "documented_not_validated"},
        {"order": 8, "stage": "render", "entrypoint_group": "config/gate12bm_script_whitelist_v1.tsv stage=render", "status": "documented_not_validated"},
        {"order": 9, "stage": "static_verify", "entrypoint_group": "workflow/run.sh verify", "status": "validated_offline"},
    ]
    write_tsv(repository / "workflow" / "stage_map.tsv", ["order", "stage", "entrypoint_group", "status"], stage_rows)

    method_rows = [
        {"evidence_block": "Discovery atlas and broad cell classes", "scripts": "scripts/run_gate3_qc.R; scripts/run_gate3b_annotation_pseudobulk.R; scripts/run_gate12a_integrated_atlas.R", "final_figure": "Figure 1; Figure S1"},
        {"evidence_block": "Myeloid state decomposition", "scripts": "scripts/run_gate12b_cell_states.R; scripts/finalize_gate12b_outputs.R", "final_figure": "Figure 2"},
        {"evidence_block": "T/NK state and transcriptional robustness", "scripts": "scripts/run_gate5b_tnk_decomposition.R; scripts/run_gate6b_confounded_robustness.R", "final_figure": "Figure 3; Figure S3"},
        {"evidence_block": "Sender-receiver hypotheses and supportive co-detection", "scripts": "scripts/run_gate12c_sender_receiver.R; scripts/run_gate12c_external_validation.R; scripts/run_gate10s_full.R", "final_figure": "Figure 4; Figure S12"},
        {"evidence_block": "External representation and Axis1", "scripts": "scripts/run_gate8a_gse266330.R; scripts/run_gate8b_oep005136.R; scripts/run_gate12v2_axis1_interpretation.R", "final_figure": "Figure 5; Figures S4-S8/S11"},
        {"evidence_block": "Spatial reconstruction, geometry, and sensitivity", "scripts": "scripts/run_gate12e_reconstruction_rctd.R; scripts/run_gate12v3_spatial_geometry.R; scripts/run_gate12ab_spatial_sensitivity_grid.R", "final_figure": "Figure 6; Figures S9-S10"},
        {"evidence_block": "Final publication rendering", "scripts": "scripts/build_gate12az_literature_standard_figure1.R; scripts/build_gate12az_literature_standard_figures23.R; scripts/build_gate12be_review_driven_figures456.R", "final_figure": "Figures 1-6"},
    ]
    write_tsv(repository / "docs" / "METHODS_TO_CODE.tsv", ["evidence_block", "scripts", "final_figure"], method_rows)

    claim_rows = [
        {"claim_scope": "Pan-cohort cellular atlas", "figure": "Figure 1", "code_group": "stage=discovery; stage=umap", "evidence_boundary": "descriptive discovery atlas"},
        {"claim_scope": "Myeloid state remodeling", "figure": "Figure 2", "code_group": "stage=discovery", "evidence_boundary": "patient-aware discovery and cohort concordance"},
        {"claim_scope": "T/NK program robustness", "figure": "Figure 3", "code_group": "stage=discovery", "evidence_boundary": "confounder-aware internal robustness"},
        {"claim_scope": "Sender-receiver hypotheses", "figure": "Figure 4", "code_group": "stage=communication; stage=supportive_projection", "evidence_boundary": "computational hypothesis plus supportive co-detection; not causal"},
        {"claim_scope": "Cross-origin representation and Axis1", "figure": "Figure 5", "code_group": "stage=external_validation", "evidence_boundary": "heterogeneous external support with preserved failures"},
        {"claim_scope": "Spatial organization", "figure": "Figure 6", "code_group": "stage=spatial", "evidence_boundary": "section-level spatial association; not patient-level causal inference"},
        {"claim_scope": "Publication rendering", "figure": "Figures 1-6 and S1-S12", "code_group": "stage=render", "evidence_boundary": "historical authoritative render chain; clean replay not yet validated"},
    ]
    write_tsv(repository / "workflow" / "claim_script_map.tsv", ["claim_scope", "figure", "code_group", "evidence_boundary"], claim_rows)

    dag_rows = [
        {"artifact": "discovery_atlas", "depends_on": "public_raw_data;annotated_discovery_sce", "availability": "partial; processed SCE not redistributed"},
        {"artifact": "lineage_states", "depends_on": "discovery_atlas", "availability": "scripts only"},
        {"artifact": "communication_hypotheses", "depends_on": "lineage_states;nichenet_priors", "availability": "priors not redistributed"},
        {"artifact": "external_axis1", "depends_on": "GSE266330;OEP005136;frozen_signatures", "availability": "restricted OEP matrices not redistributed"},
        {"artifact": "spatial_geometry", "depends_on": "GSE323357;spatial_reconstruction", "availability": "source data present; raw rerun not validated"},
        {"artifact": "final_figures", "depends_on": "discovery_atlas;lineage_states;communication_hypotheses;external_axis1;spatial_geometry", "availability": "accepted static outputs present"},
    ]
    write_tsv(repository / "workflow" / "artifact_dag.tsv", ["artifact", "depends_on", "availability"], dag_rows)

    reuse_rows = [
        {"material": "Original project code", "included": "yes; sanitized allow-list", "redistribution_status": "AUTHOR_LICENSE_DECISION_REQUIRED"},
        {"material": "Accepted figures/source data", "included": "yes; companion", "redistribution_status": "AUTHOR_CLEARANCE_REQUIRED"},
        {"material": "Raw GEO data", "included": "no", "redistribution_status": "ACCESSION_ONLY"},
        {"material": "NODE/OEP matrices", "included": "no", "redistribution_status": "PROVIDER_CONTROLLED"},
        {"material": "NicheNet priors", "included": "no", "redistribution_status": "THIRD_PARTY"},
        {"material": "Immune reference RDS", "included": "no", "redistribution_status": "THIRD_PARTY"},
        {"material": "Raw histology/source imagery", "included": "no", "redistribution_status": "THIRD_PARTY_REVIEW_REQUIRED"},
    ]
    write_tsv(repository / "LICENSES" / "DATA_REUSE_MATRIX.tsv", ["material", "included", "redistribution_status"], reuse_rows)


def build_blocker_ledger(source: Path, destination: Path) -> tuple[int, int]:
    rows = read_tsv(source)
    fields = list(rows[0]) + ["gate12bm_status", "gate12bm_resolved", "gate12bm_evidence"]
    for row in rows:
        row["gate12bm_status"] = row["gate12bl_status"]
        row["gate12bm_resolved"] = row["gate12bl_resolved"]
        row["gate12bm_evidence"] = "Inherited from the frozen Gate12BL disposition."
    new_row = {field: "" for field in fields}
    new_row.update(
        {
            "blocker_id": "BI-18",
            "category": "data availability",
            "severity": "BLOCKING",
            "finding": "Supportive communication datasets GSE225209 and GSE190772 are not fully disclosed in the active Data availability paragraph.",
            "evidence": "The active Results/Methods use GSE225209 and retain GSE190772 as non-evaluable, while Data availability lists only the five core sources.",
            "required_action": "Add both accessions and their evidence roles to the manuscript; explicitly state that GSE190772 was non-evaluable under the prespecified eligibility rule.",
            "gate12bi_resolved": "False",
            "gate12bk_status": "NOT_ASSESSED",
            "gate12bl_status": "NOT_ASSESSED",
            "gate12bl_resolved": "False",
            "gate12bm_status": "UNRESOLVED",
            "gate12bm_resolved": "False",
            "gate12bm_evidence": "Gate12BM disclosure audit; manuscript left unchanged by design.",
        }
    )
    rows.append(new_row)
    write_tsv(destination, fields, rows)
    unresolved = sum(row["gate12bm_resolved"].lower() != "true" for row in rows)
    blocking = sum(row["severity"] == "BLOCKING" and row["gate12bm_resolved"].lower() != "true" for row in rows)
    return unresolved, blocking


def assert_disclosure_gap(manuscript: Path) -> None:
    text = manuscript.read_text(encoding="utf-8")
    if "GSE225209" not in text or "GSE190772" not in text:
        raise RuntimeError("Expected supportive dataset records are absent from the active manuscript")
    match = re.search(r"(?s)^#+ Data availability\s*(.*?)(?=^#+ |\Z)", text, flags=re.MULTILINE)
    if not match:
        raise RuntimeError("Data availability section not found")
    section = match.group(1)
    if "GSE225209" in section or "GSE190772" in section:
        raise RuntimeError("Frozen BI-18 assumption changed; re-review Gate12BM rather than silently continuing")


def build(args: argparse.Namespace) -> Path:
    project = args.project_root.resolve()
    config_path = (project / args.config).resolve()
    config = read_key_value_tsv(config_path)
    output = args.output.resolve() if args.output else (project / config["output_dir"]).resolve()

    if args.smoke:
        tmp_root = Path(tempfile.gettempdir()).resolve()
        if output == tmp_root or tmp_root not in output.parents:
            raise RuntimeError("Smoke output must be inside the system temporary directory")
    elif output != (project / config["output_dir"]).resolve():
        raise RuntimeError("Formal output must match the frozen config")
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output}")

    plan = project / config["gate_plan"]
    script_allowlist_path = project / config["script_allowlist"]
    config_allowlist_path = project / config["config_allowlist"]
    builder = project / config["builder_script"]
    bl = project / config["source_gate12bl_dir"]
    receipt = bl / "GATE12BL_RECEIPT.json"
    bl_manifest = bl / "provenance" / "GATE12BL_OUTPUT_MANIFEST.tsv"
    manuscript = bl / "manuscript" / "Gate12BL_Reproducibility_Aligned_Manuscript.md"
    workbook = bl / "tables" / "Supplementary_Tables_S4_and_S6.xlsx"

    assert_hash(plan, config["expected_plan_sha256"], "Gate12BM plan")
    assert_hash(script_allowlist_path, config["expected_script_allowlist_sha256"], "script allow-list")
    assert_hash(config_allowlist_path, config["expected_config_allowlist_sha256"], "config allow-list")
    assert_hash(builder, config["expected_builder_sha256"], "Gate12BM builder")
    assert_hash(receipt, config["expected_gate12bl_receipt_sha256"], "Gate12BL receipt")
    assert_hash(bl_manifest, config["expected_gate12bl_output_manifest_sha256"], "Gate12BL output manifest")
    assert_hash(manuscript, config["expected_gate12bl_manuscript_sha256"], "Gate12BL manuscript")
    assert_hash(workbook, config["expected_supplementary_workbook_sha256"], "Gate12BL supplementary workbook")
    assert_disclosure_gap(manuscript)

    script_rows = read_tsv(script_allowlist_path)
    config_rows = read_tsv(config_allowlist_path)
    validate_allowlist(script_rows, project, "script")
    validate_allowlist(config_rows, project, "config")
    if allowlist_content_digest(script_rows, project) != config["expected_script_content_digest"]:
        raise RuntimeError("Allow-listed script content digest mismatch")
    config_digest = allowlist_content_digest(
        config_rows,
        project,
        excluded_sources={config["config_self_path"]},
    )
    if config_digest != config["expected_config_content_digest"]:
        raise RuntimeError("Allow-listed config content digest mismatch")
    upstream_manifest_rows = verify_gate12bl_manifest(bl, bl_manifest)

    stage = output.with_name(output.name + ".staging")
    if args.smoke and Path(tempfile.gettempdir()).resolve() not in stage.parents:
        raise RuntimeError("Smoke staging directory escaped the system temporary directory")
    if stage.exists():
        raise FileExistsError(f"Refusing stale staging directory: {stage}")
    stage.mkdir(parents=True)
    repository = stage / "repository"
    companion = stage / "companion"
    admin = stage / "admin"
    repository.mkdir()
    companion.mkdir()
    admin.mkdir()

    try:
        transformation_rows: list[dict[str, object]] = []
        copy_sanitized_allowlist(project, repository, script_rows, transformation_rows)
        copy_sanitized_allowlist(project, repository, config_rows, transformation_rows)
        write_text(repository / "docs" / "GATE12BM_PLAN.md", plan.read_text(encoding="utf-8"))
        build_repository_documents(repository)
        build_ledgers(repository)

        copy_tree_exact(bl / "figures", companion / "reference" / "figures")
        copy_tree_exact(bl / "source_data", companion / "data" / "source_data")
        copy_tree_exact(bl / "tables", companion / "reference" / "tables")
        copy_tree_exact(bl / "manuscript", companion / "reference" / "manuscript")
        write_text(
            companion / "README.md",
            """# Static evidence companion candidate

This directory contains exact, hash-verified copies of the accepted Gate12BL
figures, source-data files, supplementary tables, manuscript, and legends. It is
intended as a local Zenodo-style companion candidate. Inclusion here does not
constitute author approval to redistribute; license and third-party reuse review
remain pending. Restricted raw matrices and third-party reference objects are not
included.
""",
        )
        original_receipt = json.loads(receipt.read_text(encoding="utf-8"))
        sanitized_receipt = {
            "source_gate": "Gate12BL run_v1",
            "source_receipt_sha256": sha256_file(receipt),
            "source_status": original_receipt.get("status", "UNKNOWN"),
            "source_verdicts": original_receipt.get("verdicts", {}),
            "active_manuscript_sha256": sha256_file(manuscript),
            "supplementary_workbook_sha256": sha256_file(workbook),
            "figure_files": 42,
            "source_data_files": 85,
            "scope": "sanitized upstream identity only; original path-bearing receipt excluded",
        }
        write_text(companion / "reference" / "receipts" / "GATE12BL_PUBLIC_SOURCE_RECEIPT.json", json.dumps(sanitized_receipt, indent=2, sort_keys=True))

        figure_count = len([p for p in (companion / "reference" / "figures").rglob("*") if p.is_file()])
        source_data_count = len([p for p in (companion / "data" / "source_data").rglob("*") if p.is_file()])
        table_files = sorted(p.name for p in (companion / "reference" / "tables").glob("*") if p.is_file())
        manuscript_files = sorted(p.name for p in (companion / "reference" / "manuscript").glob("*") if p.is_file())
        if figure_count != int(config["expected_figures"]):
            raise RuntimeError(f"Expected {config['expected_figures']} figure files, observed {figure_count}")
        if source_data_count != int(config["expected_source_data_files"]):
            raise RuntimeError(f"Expected {config['expected_source_data_files']} source-data files, observed {source_data_count}")
        expected_tables = sorted(config["expected_table_files"].split(";"))
        expected_manuscripts = sorted(config["expected_manuscript_files"].split(";"))
        if table_files != expected_tables:
            raise RuntimeError(f"Table payload mismatch: expected {expected_tables}, observed {table_files}")
        if manuscript_files != expected_manuscripts:
            raise RuntimeError(f"Manuscript payload mismatch: expected {expected_manuscripts}, observed {manuscript_files}")

        script_transform_rows = [row for row in transformation_rows if row["candidate_path"].startswith("repository/scripts/")]
        config_transform_rows = [row for row in transformation_rows if row["candidate_path"].startswith("repository/config/")]
        write_tsv(
            repository / "manifests" / "scripts.tsv",
            ["source_path", "candidate_path", "stage", "declared_status", "transformations", "source_sha256", "candidate_sha256", "byte_identical"],
            script_transform_rows,
        )
        write_tsv(
            repository / "manifests" / "configs.tsv",
            ["source_path", "candidate_path", "stage", "declared_status", "transformations", "source_sha256", "candidate_sha256", "byte_identical"],
            config_transform_rows,
        )
        data_role_rows = read_tsv(repository / "data" / "accessions.tsv")
        write_tsv(
            repository / "manifests" / "data_roles.tsv",
            ["accession_or_id", "provider", "role", "access", "url", "redistributed"],
            data_role_rows,
        )
        expected_output_rows = [
            {"artifact": "accepted_figure_files", "expected": figure_count, "validation": "hash_verified_static_payload"},
            {"artifact": "accepted_source_data_files", "expected": source_data_count, "validation": "hash_verified_static_payload"},
            {"artifact": "accepted_table_files", "expected": len(table_files), "validation": "hash_verified_static_payload"},
            {"artifact": "active_manuscript_and_legends", "expected": len(manuscript_files), "validation": "hash_verified_static_payload"},
            {"artifact": "offline_verify", "expected": "PASS", "validation": "validated_in_gate12bm"},
            {"artifact": "figure_replay", "expected": "NOT_VALIDATED", "validation": "future_gate_required"},
            {"artifact": "full_raw_reconstruction", "expected": "NOT_VALIDATED", "validation": "future_gate_required"},
        ]
        write_tsv(repository / "manifests" / "expected_outputs.tsv", ["artifact", "expected", "validation"], expected_output_rows)

        syntax_rows = syntax_audit(repository, config.get("rscript", "Rscript"))
        security_rows = security_audit(repository, companion)
        write_tsv(admin / "GATE12BM_TRANSFORMATION_LEDGER.tsv", ["source_path", "candidate_path", "stage", "declared_status", "transformations", "source_sha256", "candidate_sha256", "byte_identical"], transformation_rows)
        write_tsv(admin / "GATE12BM_CODE_SYNTAX_AUDIT.tsv", ["relative_path", "language", "status", "detail"], syntax_rows)
        write_tsv(admin / "GATE12BM_SECURITY_PRIVACY_AUDIT.tsv", ["scope", "path", "rule", "status", "detail"], security_rows)
        unresolved, blocking = build_blocker_ledger(bl / "admin" / "GATE12BL_BLOCKER_DISPOSITION.tsv", admin / "GATE12BM_BLOCKER_DISPOSITION.tsv")

        repository_manifest = repository / "manifests" / "repository_files.tsv"
        companion_manifest = companion / "manifests" / "companion_files.tsv"
        write_tsv(repository_manifest, ["relative_path", "bytes", "sha256"], manifest_rows(repository, {repository_manifest}))
        write_tsv(companion_manifest, ["relative_path", "bytes", "sha256"], manifest_rows(companion, {companion_manifest}))
        verify_manifest(repository, repository_manifest)
        verify_manifest(companion, companion_manifest)

        verification = subprocess.run(
            [str(repository / "workflow" / "run.sh"), "verify"],
            cwd=repository,
            text=True,
            capture_output=True,
            check=False,
        )
        if verification.returncode != 0:
            raise RuntimeError(f"Portable verifier failed: {verification.stderr or verification.stdout}")
        write_text(admin / "GATE12BM_OFFLINE_VERIFY.txt", verification.stdout + verification.stderr)

        acceptance_rows = [
            {"dimension": "upstream_identity", "status": "PASS", "evidence": "Frozen Gate12BL receipt/manuscript/workbook hashes match"},
            {"dimension": "upstream_selected_payload_closure", "status": "PASS", "evidence": f"{upstream_manifest_rows} Gate12BL manifest rows verified; selected publication directories are closed against that manifest"},
            {"dimension": "package_structure", "status": "PASS", "evidence": "Repository, companion and admin payloads built"},
            {"dimension": "manifest_integrity", "status": "PASS", "evidence": "Repository and companion manifests independently rehashed"},
            {"dimension": "static_publication_payload", "status": "PASS", "evidence": f"{figure_count} figure files and {source_data_count} source-data files"},
            {"dimension": "code_syntax", "status": "PASS", "evidence": f"{len(syntax_rows)} public scripts parsed"},
            {"dimension": "privacy_security", "status": "PASS", "evidence": "No finding under the defined automated repository/companion rules"},
            {"dimension": "repository_size_type_policy", "status": "PASS", "evidence": "No disallowed binary/raw object, file >5 MiB, or symlink"},
            {"dimension": "companion_boundary", "status": "PASS", "evidence": "Large/static evidence isolated from code repository"},
            {"dimension": "accession_ledger", "status": "PASS", "evidence": "Core, supportive, non-evaluable, and provider-controlled roles documented"},
            {"dimension": "offline_verification", "status": "PASS", "evidence": "workflow/run.sh verify returned zero"},
            {"dimension": "archive_roundtrip", "status": "PASS", "evidence": "Both archives extracted to a clean temporary directory and bundled verify returned zero"},
            {"dimension": "figure_replay", "status": "NOT_VALIDATED", "evidence": "No clean rerender performed in Gate12BM"},
            {"dimension": "clean_environment_smoke", "status": "NOT_VALIDATED", "evidence": "No validated lockfile/container"},
            {"dimension": "full_raw_reconstruction", "status": "NOT_VALIDATED", "evidence": "Restricted and third-party inputs are not bundled"},
            {"dimension": "redistribution_clearance", "status": "PENDING_AUTHOR_APPROVAL", "evidence": "Author/third-party review required"},
            {"dimension": "repository_license", "status": "NOT_SELECTED", "evidence": "No license fabricated"},
            {"dimension": "persistent_identifier", "status": "NOT_ASSIGNED", "evidence": "No remote upload or DOI minting performed"},
            {"dimension": "manuscript_dataset_disclosure", "status": "FAIL", "evidence": "BI-18 unresolved for GSE225209/GSE190772"},
            {"dimension": "submission_readiness", "status": "NOT_READY", "evidence": f"{blocking} unresolved blocking items"},
        ]
        write_tsv(admin / "GATE12BM_ACCEPTANCE_MATRIX.tsv", ["dimension", "status", "evidence"], acceptance_rows)

        report = f"""# Gate12BM public repository candidate report

- Verdict: `{EXPECTED_VERDICT}`
- Formal execution: `{str(not args.smoke).upper()}`
- Accepted figure files: {figure_count}
- Accepted source-data files: {source_data_count}
- Sanitized allow-listed files: {len(transformation_rows)}
- Parsed public scripts: {len(syntax_rows)}
- Unresolved blockers: {unresolved}
- Unresolved blocking items: {blocking}
- Static source-data verification: `PASS`
- Figure replay: `NOT_VALIDATED`
- Clean-environment smoke: `NOT_VALIDATED`
- Full raw reconstruction: `NOT_VALIDATED`
- Redistribution clearance: `PENDING_AUTHOR_APPROVAL`
- Repository license: `NOT_SELECTED`
- Persistent identifier: `NOT_ASSIGNED`
- Submission readiness: `NOT_READY`

Manifest convention: repository and companion manifests exclude themselves; the
top-level output manifest excludes itself and `GATE12BM_RECEIPT.json` to avoid
self-hash cycles. The receipt records the output-manifest hash, and the final
receipt hash must be reported outside the package after execution.

Gate12BM created a sanitized local candidate only. It did not modify the active
manuscript, upload files, create a remote repository, choose a license, or mint a
DOI. BI-18 records the required disclosure repair for GSE225209 and GSE190772.
"""
        write_text(stage / "GATE12BM_REPORT.md", report)

        archives = stage / "archives"
        github_archive = archives / "bone-metastasis-github-candidate.tar.gz"
        companion_archive = archives / "bone-metastasis-zenodo-companion-candidate.tar.gz"
        # The two archives intentionally extract to sibling `repository/` and
        # `companion/` directories so the offline verifier works unchanged.
        make_deterministic_tar_gz(repository, github_archive, "repository")
        make_deterministic_tar_gz(companion, companion_archive, "companion")
        archive_roundtrip_log = archive_roundtrip_verify(github_archive, companion_archive)
        write_text(admin / "GATE12BM_ARCHIVE_ROUNDTRIP_VERIFY.txt", archive_roundtrip_log)

        output_manifest = admin / "GATE12BM_OUTPUT_MANIFEST.tsv"
        receipt_path = stage / "GATE12BM_RECEIPT.json"
        write_tsv(
            output_manifest,
            ["relative_path", "bytes", "sha256"],
            manifest_rows(stage, {output_manifest, receipt_path}),
        )
        verify_manifest(stage, output_manifest)
        actual_output_files = {
            path.relative_to(stage).as_posix()
            for path in stage.rglob("*")
            if path.is_file() and path not in {output_manifest, receipt_path}
        }
        listed_output_files = {row["relative_path"] for row in read_tsv(output_manifest)}
        if actual_output_files != listed_output_files:
            raise RuntimeError("Gate12BM output manifest is not closed over the staged payload")
        receipt_payload = {
            "schema_version": SCHEMA_VERSION,
            "generated_at_utc": utc_now(),
            "formal_execution": not args.smoke,
            "verdict": EXPECTED_VERDICT,
            "submission_readiness": "NOT_READY",
            "source_gate": "Gate12BL run_v1",
            "source_gate_receipt_sha256": sha256_file(receipt),
            "source_gate_output_manifest_sha256": sha256_file(bl_manifest),
            "source_gate_manifest_rows_verified": upstream_manifest_rows,
            "active_manuscript_sha256": sha256_file(manuscript),
            "figures": figure_count,
            "source_data_files": source_data_count,
            "allowlisted_scripts": len(script_rows),
            "allowlisted_configs": len(config_rows),
            "syntax_checked_scripts": len(syntax_rows),
            "unresolved_blockers": unresolved,
            "unresolved_blocking_items": blocking,
            "new_blocker": "BI-18",
            "figure_replay": "NOT_VALIDATED",
            "clean_environment_smoke": "NOT_VALIDATED",
            "full_raw_reconstruction": "NOT_VALIDATED",
            "redistribution_clearance": "PENDING_AUTHOR_APPROVAL",
            "repository_license": "NOT_SELECTED",
            "persistent_identifier": "NOT_ASSIGNED",
            "repository_manifest_sha256": sha256_file(repository_manifest),
            "companion_manifest_sha256": sha256_file(companion_manifest),
            "output_manifest_sha256": sha256_file(output_manifest),
            "output_manifest_exclusions": [
                "admin/GATE12BM_OUTPUT_MANIFEST.tsv",
                "GATE12BM_RECEIPT.json",
            ],
            "github_candidate_archive_sha256": sha256_file(github_archive),
            "companion_candidate_archive_sha256": sha256_file(companion_archive),
            "external_mutations": [],
        }
        write_text(receipt_path, json.dumps(receipt_payload, indent=2, sort_keys=True))

        stage.rename(output)
        return output
    except Exception:
        if stage.exists():
            shutil.rmtree(stage)
        raise


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_project = Path(__file__).resolve().parents[1]
    parser.add_argument("--project-root", type=Path, default=default_project)
    parser.add_argument("--config", type=Path, default=Path("config/gate12bm_public_repository_candidate_v1.tsv"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output = build(args)
    print(f"Gate12BM output: {output}")
    print(f"Verdict: {EXPECTED_VERDICT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
