#!/usr/bin/env python3

"""Classify runtime packages that may require corresponding-source review."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def classify(license_text: str) -> tuple[str, str]:
    upper = license_text.upper()
    if "GPL" not in upper:
        return "notice-only-review", "No GPL/LGPL token in package metadata"
    if " OR " in upper:
        return "alternative-license-review", "Metadata declares a license alternative"
    if "EXCEPTION" in upper:
        return "exception-review", "Metadata declares a linking/runtime exception"
    if "LGPL" in upper:
        return "lgpl-source-review", "LGPL-family binary in runtime image"
    return "gpl-source-review", "GPL-family binary in runtime image"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with args.inventory.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    required = {"name", "version", "build", "license", "binary_url", "sha256"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("Inventory has no rows or required columns")

    fields = [
        "name",
        "version",
        "build",
        "license",
        "review_class",
        "collect_source",
        "review_note",
        "binary_url",
        "binary_sha256",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda item: item["name"]):
            review_class, note = classify(row["license"])
            if review_class == "notice-only-review":
                continue
            writer.writerow({
                "name": row["name"],
                "version": row["version"],
                "build": row["build"],
                "license": row["license"],
                "review_class": review_class,
                "collect_source": "REVIEW",
                "review_note": note,
                "binary_url": row["binary_url"],
                "binary_sha256": row["sha256"],
            })


if __name__ == "__main__":
    main()
