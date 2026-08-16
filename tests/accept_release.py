#!/usr/bin/env python3
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
        if ".git" in path.relative_to(root).parts:
            continue
        if path.is_symlink():
            failures.append(f"symlink:{path}")
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path != manifest and ".git" not in path.relative_to(root).parts
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
