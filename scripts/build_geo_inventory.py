#!/usr/bin/env python3
"""Build a reproducible GEO sample and supplementary-file inventory for Gate 2."""

from __future__ import annotations

import csv
import gzip
import html
import re
import sys
import time
from pathlib import Path

import requests


ACCESSIONS = [
    "GSE266330",
    "GSE143791",
    "GSE202813",
    "GSE225209",
    "GSE225208",
    "GSE190772",
    "GSE39494",
    "GSE14776",
]

USER_AGENT = "gate2-multicohort-inventory/1.0 (public GEO metadata audit)"


def geo_prefix(accession: str) -> str:
    return accession[:-3] + "nnn"


def request(session: requests.Session, url: str) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, 5):
        try:
            response = session.get(url, timeout=(20, 120))
            response.raise_for_status()
            time.sleep(0.34)
            return response
        except Exception as exc:
            last_error = exc
            if attempt < 4:
                time.sleep(2**attempt)
    raise RuntimeError(f"Failed to retrieve {url}: {last_error}")


def download_soft(session: requests.Session, accession: str, raw_dir: Path) -> Path:
    url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{geo_prefix(accession)}/{accession}/soft/{accession}_family.soft.gz"
    )
    destination = raw_dir / f"{accession}_family.soft.gz"
    response = request(session, url)
    temporary = destination.with_suffix(destination.suffix + ".part")
    temporary.write_bytes(response.content)
    if len(response.content) == 0:
        raise RuntimeError(f"Empty SOFT response for {accession}")
    temporary.replace(destination)
    return destination


def parse_soft(path: Path, accession: str) -> tuple[dict[str, str], list[dict[str, str]]]:
    series: dict[str, str] = {"accession": accession}
    samples: list[dict[str, str]] = []
    current: dict[str, object] | None = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SERIES ="):
                continue
            if line.startswith("^SAMPLE ="):
                if current is not None:
                    samples.append(flatten_sample(current, accession))
                current = {
                    "gsm": line.split("=", 1)[1].strip(),
                    "characteristics": [],
                    "supplementary_files": [],
                    "relations": [],
                }
                continue
            if current is None:
                if line.startswith("!Series_title ="):
                    series["title"] = line.split("=", 1)[1].strip()
                elif line.startswith("!Series_summary =") and "summary" not in series:
                    series["summary"] = line.split("=", 1)[1].strip()
                elif line.startswith("!Series_overall_design =") and "overall_design" not in series:
                    series["overall_design"] = line.split("=", 1)[1].strip()
                elif line.startswith("!Series_pubmed_id ="):
                    series["pubmed_id"] = line.split("=", 1)[1].strip()
                continue

            if line.startswith("!Sample_title ="):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 ="):
                current["source_name"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_organism_ch1 ="):
                current["organism"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 ="):
                current["characteristics"].append(line.split("=", 1)[1].strip())
            elif line.startswith("!Sample_supplementary_file"):
                current["supplementary_files"].append(line.split("=", 1)[1].strip())
            elif line.startswith("!Sample_relation ="):
                current["relations"].append(line.split("=", 1)[1].strip())
        if current is not None:
            samples.append(flatten_sample(current, accession))
    return series, samples


def flatten_sample(sample: dict[str, object], accession: str) -> dict[str, str]:
    return {
        "accession": accession,
        "gsm": str(sample.get("gsm", "")),
        "title": str(sample.get("title", "")),
        "source_name": str(sample.get("source_name", "")),
        "organism": str(sample.get("organism", "")),
        "characteristics": " | ".join(sample.get("characteristics", [])),
        "supplementary_files": " | ".join(sample.get("supplementary_files", [])),
        "relations": " | ".join(sample.get("relations", [])),
    }


def parse_filelist(text: str, accession: str, base_url: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        name = row.get("Name", "")
        if not name:
            continue
        rows.append(
            {
                "accession": accession,
                "name": name,
                "bytes": row.get("Size", ""),
                "type": row.get("Type", ""),
                "url": base_url + name,
                "inventory_basis": "GEO filelist.txt",
            }
        )
    return rows


def parse_directory(text: str, accession: str, base_url: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for match in re.finditer(r'<a href="([^"]+)">', text):
        name = html.unescape(match.group(1))
        if name.startswith(("http", "/", "?")) or name in {"../", "filelist.txt"}:
            continue
        rows.append(
            {
                "accession": accession,
                "name": name,
                "bytes": "",
                "type": Path(name).suffix.lstrip(".").upper(),
                "url": base_url + name,
                "inventory_basis": "GEO supplementary directory",
            }
        )
    return rows


def supplementary_inventory(
    session: requests.Session, accession: str
) -> list[dict[str, str]]:
    base_url = (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/"
        f"{geo_prefix(accession)}/{accession}/suppl/"
    )
    directory = request(session, base_url)
    if "filelist.txt" in directory.text:
        filelist = request(session, base_url + "filelist.txt")
        return parse_filelist(filelist.text, accession, base_url)
    return parse_directory(directory.text, accession, base_url)


def write_tsv(rows: list[dict[str, str]], path: Path, fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    output = Path(sys.argv[1] if len(sys.argv) > 1 else "metadata").resolve()
    raw_dir = output / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    series_rows: list[dict[str, str]] = []
    sample_rows: list[dict[str, str]] = []
    file_rows: list[dict[str, str]] = []
    for accession in ACCESSIONS:
        print(f"INVENTORY\t{accession}", flush=True)
        soft = download_soft(session, accession, raw_dir)
        series, samples = parse_soft(soft, accession)
        series_rows.append(series)
        sample_rows.extend(samples)
        file_rows.extend(supplementary_inventory(session, accession))

    write_tsv(
        series_rows,
        output / "geo_series.tsv",
        ["accession", "title", "pubmed_id", "summary", "overall_design"],
    )
    write_tsv(
        sample_rows,
        output / "geo_samples.tsv",
        [
            "accession",
            "gsm",
            "title",
            "source_name",
            "organism",
            "characteristics",
            "supplementary_files",
            "relations",
        ],
    )
    write_tsv(
        file_rows,
        output / "geo_files.tsv",
        ["accession", "name", "bytes", "type", "url", "inventory_basis"],
    )
    print(
        f"COMPLETE\tseries={len(series_rows)}\tsamples={len(sample_rows)}\tfiles={len(file_rows)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
