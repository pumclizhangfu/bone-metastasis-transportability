#!/usr/bin/env python3
"""Verify that every manifest file has the expected size and checksum receipt."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    cfg = parser.parse_args()
    with cfg.manifest.open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    failures = 0
    for item in manifest:
        path = cfg.data / item["accession"] / item["filename"]
        checksum_path = Path(str(path) + ".sha256")
        exists = path.is_file()
        observed = path.stat().st_size if exists else -1
        expected = int(item["bytes"])
        checksum = ""
        if checksum_path.is_file():
            checksum = checksum_path.read_text(encoding="utf-8").split()[0]
        status = "PASS" if exists and observed == expected and len(checksum) == 64 else "FAIL"
        failures += status == "FAIL"
        rows.append(
            {
                "accession": item["accession"],
                "filename": item["filename"],
                "expected_bytes": expected,
                "observed_bytes": observed,
                "sha256": checksum,
                "integrity_status": status,
                "role": item["role"],
            }
        )
    cfg.output.parent.mkdir(parents=True, exist_ok=True)
    with cfg.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"COMPLETE\tfiles={len(rows)}\tfailures={failures}\toutput={cfg.output}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
