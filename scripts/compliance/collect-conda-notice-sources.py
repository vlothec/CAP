#!/usr/bin/env python3

"""Collect reviewed supplemental Conda notice source files."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import urllib.request


def digest(path: Path, algorithm: str) -> str:
    if algorithm not in {"md5", "sha256"}:
        raise ValueError(f"Unsupported digest type: {algorithm}")
    checksum = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reviewed_sources", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")
    with args.reviewed_sources.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "package", "version", "build", "source_url", "digest_type",
        "expected_digest", "source_role",
    }
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("Reviewed source CSV has no rows or required columns")

    args.output.mkdir(parents=True)
    collected: list[dict[str, str]] = []
    counters: dict[tuple[str, str, str], int] = {}
    for row in rows:
        key = (row["package"], row["version"], row["build"])
        counters[key] = counters.get(key, 0) + 1
        filename = row["source_url"].rsplit("/", 1)[-1] or "source"
        directory = args.output / "upstream-source" / "-".join(key)
        path = directory / f"{counters[key]:02d}-{filename}"
        path.parent.mkdir(parents=True, exist_ok=True)
        request = urllib.request.Request(
            row["source_url"], headers={"User-Agent": "CAP-compliance/1"}
        )
        with urllib.request.urlopen(request) as response:
            path.write_bytes(response.read())
        actual = digest(path, row["digest_type"])
        if actual.lower() != row["expected_digest"].lower():
            raise SystemExit(f"Source checksum mismatch: {row['source_url']}")
        collected.append({
            **row,
            "collected_path": path.relative_to(args.output).as_posix(),
            "sha256": digest(path, "sha256"),
        })

    fields = list(collected[0])
    with (args.output / "source-to-binary.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(collected)
    with (args.output / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in sorted(p for p in args.output.rglob("*") if p.is_file()
                           and p.name != "SHA256SUMS"):
            handle.write(
                f"{digest(path, 'sha256')}  {path.relative_to(args.output)}\n"
            )
    print(f"Collected {len(collected)} supplemental source files")


if __name__ == "__main__":
    main()
