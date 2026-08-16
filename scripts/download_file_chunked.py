#!/usr/bin/env python3
"""Download one public file with resumable, verified HTTP-range chunks."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import tarfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests


LOCAL = threading.local()


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--chunk-mib", type=int, default=8)
    parser.add_argument("--attempts", type=int, default=8)
    return parser.parse_args()


def atomic_json(path: Path, record: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(record, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_range(path: Path, start: int, response: requests.Response, expected: int) -> str:
    digest = hashlib.sha256()
    written = 0
    with path.open("r+b", buffering=0) as handle:
        handle.seek(start)
        for block in response.iter_content(1024 * 1024):
            if not block:
                continue
            if written + len(block) > expected:
                raise RuntimeError("server returned more bytes than requested")
            handle.write(block)
            digest.update(block)
            written += len(block)
        handle.flush()
        os.fsync(handle.fileno())
    if written != expected:
        raise EOFError(f"truncated range: {written}/{expected}")
    return digest.hexdigest()


def fetch(
    url: str,
    output: Path,
    state_dir: Path,
    index: int,
    start: int,
    end: int,
    total: int,
    attempts: int,
) -> tuple[int, int]:
    marker = state_dir / f"chunk_{index:06d}.json"
    expected = end - start + 1
    if marker.exists():
        return index, expected
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            session = getattr(LOCAL, "session", None)
            if session is None:
                session = requests.Session()
                session.headers.update({"User-Agent": "bone-metastasis-gate2/1.0"})
                LOCAL.session = session
            with session.get(
                url,
                headers={"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity"},
                stream=True,
                timeout=(30, 180),
            ) as response:
                if response.status_code != 206:
                    raise RuntimeError(f"expected HTTP 206, received {response.status_code}")
                wanted = f"bytes {start}-{end}/{total}"
                if response.headers.get("Content-Range") != wanted:
                    raise RuntimeError(
                        f"Content-Range mismatch: {response.headers.get('Content-Range')!r} != {wanted!r}"
                    )
                digest = write_range(output, start, response, expected)
            atomic_json(
                marker,
                {"index": index, "start": start, "end": end, "bytes": expected, "sha256": digest},
            )
            return index, expected
        except Exception as exc:
            last_error = exc
            print(f"RETRY\tchunk={index}\tattempt={attempt}\t{exc}", file=sys.stderr, flush=True)
            if attempt < attempts:
                time.sleep(min(30, 2**attempt))
    raise RuntimeError(f"chunk {index} failed: {last_error}")


def integrity_check(path: Path) -> str:
    if path.name.endswith(".tar"):
        with tarfile.open(path, "r") as archive:
            members = archive.getmembers()
        if not members:
            raise RuntimeError("tar archive contains no members")
        return f"tar_members={len(members)}"
    if path.name.endswith(".gz"):
        uncompressed = 0
        with gzip.open(path, "rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                uncompressed += len(block)
        if uncompressed == 0:
            raise RuntimeError("gzip file decompresses to zero bytes")
        return f"gzip_uncompressed={uncompressed}"
    return "size_only"


def main() -> int:
    cfg = args()
    if not 1 <= cfg.workers <= 32:
        raise ValueError("workers must be 1..32")
    output = cfg.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    state_dir = Path(str(output) + ".chunks")
    state_dir.mkdir(parents=True, exist_ok=True)
    state_file = state_dir / "state.json"
    chunk_size = cfg.chunk_mib * 1024 * 1024
    state = {"url": cfg.url, "size": cfg.size, "chunk_size": chunk_size}
    if state_file.exists() and json.loads(state_file.read_text()) != state:
        raise RuntimeError("existing state does not match requested download")
    if not state_file.exists():
        atomic_json(state_file, state)

    if output.exists() and output.stat().st_size not in {0, cfg.size}:
        raise RuntimeError(f"unexpected existing size: {output.stat().st_size}")
    with output.open("a+b") as handle:
        handle.truncate(cfg.size)

    ranges = []
    for index, start in enumerate(range(0, cfg.size, chunk_size)):
        end = min(cfg.size - 1, start + chunk_size - 1)
        if not (state_dir / f"chunk_{index:06d}.json").exists():
            ranges.append((index, start, end))
    total_chunks = (cfg.size + chunk_size - 1) // chunk_size
    completed = total_chunks - len(ranges)
    print(
        f"FILE\t{output.name}\tbytes={cfg.size}\tchunks={total_chunks}\t"
        f"existing={completed}\tworkers={cfg.workers}",
        flush=True,
    )
    with ThreadPoolExecutor(max_workers=cfg.workers) as pool:
        futures = [
            pool.submit(
                fetch,
                cfg.url,
                output,
                state_dir,
                index,
                start,
                end,
                cfg.size,
                cfg.attempts,
            )
            for index, start, end in ranges
        ]
        for future in as_completed(futures):
            index, _ = future.result()
            completed += 1
            print(f"CHUNK_OK\t{index}\tprogress={completed}/{total_chunks}", flush=True)

    if output.stat().st_size != cfg.size:
        raise RuntimeError("final file size mismatch")
    detail = integrity_check(output)
    digest = sha256sum(output)
    Path(str(output) + ".sha256").write_text(f"{digest}  {output.name}\n", encoding="utf-8")
    print(f"COMPLETE\t{output}\tsha256={digest}\t{detail}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
