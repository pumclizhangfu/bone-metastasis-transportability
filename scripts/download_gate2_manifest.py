#!/usr/bin/env python3
"""Download the Gate 2 manifest sequentially with eight resumable ranges per file."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=ROOT / "metadata/gate2_download_manifest.tsv")
    parser.add_argument("--output", type=Path, default=ROOT / "data/raw")
    parser.add_argument("--workers", type=int, default=8)
    cfg = parser.parse_args()
    cfg.output.mkdir(parents=True, exist_ok=True)
    with cfg.manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    total = sum(int(row["bytes"]) for row in rows)
    print(f"MANIFEST\tfiles={len(rows)}\tbytes={total}\tworkers={cfg.workers}", flush=True)
    for index, row in enumerate(rows, 1):
        destination = cfg.output / row["accession"] / row["filename"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        print(f"START\t{index}/{len(rows)}\t{row['accession']}\t{row['filename']}", flush=True)
        command = [
            sys.executable,
            str(ROOT / "scripts/download_file_chunked.py"),
            "--url",
            row["url"],
            "--output",
            str(destination),
            "--size",
            row["bytes"],
            "--workers",
            str(cfg.workers),
        ]
        subprocess.run(command, check=True)
    print("MANIFEST_COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
