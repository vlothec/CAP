#!/usr/bin/env python3

from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import tempfile

script = Path(__file__).with_name("merge-conda-source-corpora.py")
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    primary = root / "primary"
    supplemental = root / "supplemental"
    output = root / "output"
    (primary / "upstream-source").mkdir(parents=True)
    (supplemental / "upstream-source").mkdir(parents=True)
    (primary / "upstream-source/a.tar.gz").write_text("a")
    (supplemental / "upstream-source/b.tar.gz").write_text("b")
    fields = [
        "package", "version", "build", "binary_url", "binary_sha256",
        "recipe_repository", "recipe_commit", "source_url",
        "recipe_hash_algorithm", "recipe_hash_value", "source_sha256",
        "target_directory", "collected_path", "review_status",
    ]
    with (primary / "source-to-binary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow({
            "package": "a", "version": "1", "build": "0",
            "collected_path": "upstream-source/a.tar.gz",
        })
    supplemental_fields = [
        "package", "version", "build", "source_url", "digest_type",
        "expected_digest", "source_role", "collected_path", "sha256",
    ]
    with (supplemental / "source-to-binary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=supplemental_fields)
        writer.writeheader()
        writer.writerow({
            "package": "b", "version": "2", "build": "0",
            "digest_type": "sha256", "source_role": "supplemental",
            "collected_path": "upstream-source/b.tar.gz",
        })
    subprocess.run([str(script), str(primary), str(supplemental), str(output)], check=True)
    with (output / "source-to-binary.csv").open(newline="") as handle:
        result = list(csv.DictReader(handle))
    assert result[0]["collected_path"] == "primary/upstream-source/a.tar.gz"
    assert result[1]["collected_path"] == "supplemental/upstream-source/b.tar.gz"
    for row in result:
        assert (output / row["collected_path"]).is_file()
    assert not primary.exists()
    assert not supplemental.exists()
print("Conda source-corpus merge test passed.")
