#!/usr/bin/env python3

"""Merge primary and supplemental verified Conda source corpora."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import shutil


def read_rows(
    path: Path, required: set[str]
) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
        if not rows or not required.issubset(fields):
            raise SystemExit(f"Source mapping lacks required rows or columns: {path}")
        return fields, rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("primary", type=Path)
    parser.add_argument("supplemental", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")

    common = {"package", "version", "build", "collected_path"}
    fields, primary_rows = read_rows(
        args.primary / "source-to-binary.csv", common
    )
    _, supplemental_rows = read_rows(
        args.supplemental / "source-to-binary.csv",
        common | {
            "source_url", "digest_type", "expected_digest", "source_role",
            "sha256",
        },
    )
    args.output.mkdir(parents=True)
    shutil.move(args.primary, args.output / "primary")
    shutil.move(args.supplemental, args.output / "supplemental")

    rows: list[dict[str, str]] = []
    for row in primary_rows:
        source = args.output / "primary" / row["collected_path"]
        if not source.is_file():
            raise SystemExit(f"Primary source file is missing: {source}")
        row["collected_path"] = f"primary/{row['collected_path']}"
        rows.append(row)
    for source in supplemental_rows:
        source_path = args.output / "supplemental" / source["collected_path"]
        if not source_path.is_file():
            raise SystemExit(f"Supplemental source file is missing: {source_path}")
        rows.append({
            "package": source["package"],
            "version": source["version"],
            "build": source["build"],
            "source_url": source["source_url"],
            "recipe_hash_algorithm": source["digest_type"],
            "recipe_hash_value": source["expected_digest"],
            "source_sha256": source["sha256"],
            "collected_path": f"supplemental/{source['collected_path']}",
            "review_status": source["source_role"],
        })

    output = args.output / "source-to-binary.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Merged {len(primary_rows)} primary and {len(supplemental_rows)} "
        "supplemental Conda source mappings"
    )


if __name__ == "__main__":
    main()
