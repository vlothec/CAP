#!/usr/bin/env python3

"""Test conservative source-review policy application and drift detection."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import subprocess
import tempfile

PROJECT = Path(__file__).resolve().parents[2]
SCRIPT = PROJECT / "scripts" / "compliance" / "apply-source-review-policy.py"
spec = importlib.util.spec_from_file_location("source_policy", SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

fields = [
    "name", "version", "build", "license", "review_class",
    "collect_source", "review_note", "binary_url", "binary_sha256",
]
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    source = root / "review.csv"
    output = root / "reviewed.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for review_class, names in module.EXPECTED.items():
            for name in sorted(names):
                writer.writerow({
                    "name": name, "version": "1", "build": "0",
                    "license": "test", "review_class": review_class,
                    "collect_source": "REVIEW", "review_note": "test",
                    "binary_url": "https://example.invalid/package",
                    "binary_sha256": "a" * 64,
                })
    subprocess.run([str(SCRIPT), str(source), str(output)], check=True)
    with output.open(newline="", encoding="utf-8") as handle:
        reviewed = list(csv.DictReader(handle))
    assert len(reviewed) == 53
    assert all(row["collect_source"] == "yes" for row in reviewed)
    assert all("preserve exact" in row["review_note"] for row in reviewed)

    lines = source.read_text().splitlines()
    source.write_text("\n".join(lines[:-1]) + "\n")
    failure = subprocess.run(
        [str(SCRIPT), str(source), str(root / "should-fail.csv")],
        text=True, capture_output=True,
    )
    assert failure.returncode != 0
    assert "Candidate set changed" in failure.stderr

print("Source-review policy test passed.")
