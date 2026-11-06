#!/usr/bin/env python3

"""Collect embedded license evidence from exact Conda binary artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
import zipfile

LICENSE_NAME = re.compile(
    r"^(copying|copyright|licen[cs]e|notice|authors?)([._-].*)?$",
    re.IGNORECASE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "CAP-compliance/1"})
    with urllib.request.urlopen(request) as response, destination.open("wb") as output:
        shutil.copyfileobj(response, output)


def wanted(member: tarfile.TarInfo) -> bool:
    path = Path(member.name)
    if path.is_absolute() or ".." in path.parts or not member.isfile():
        return False
    return (
        member.name in {"info/about.json", "info/index.json"}
        or member.name.startswith("info/recipe/")
        or LICENSE_NAME.match(path.name) is not None
    )


def extract_tar(archive: tarfile.TarFile, destination: Path) -> None:
    archive.extractall(
        destination,
        members=[member for member in archive if wanted(member)],
        filter="data",
    )


def extract_evidence(artifact: Path, destination: Path) -> None:
    destination.mkdir(parents=True)
    if artifact.name.endswith(".tar.bz2"):
        with tarfile.open(artifact, "r:bz2") as archive:
            extract_tar(archive, destination)
        return
    if not artifact.name.endswith(".conda"):
        raise ValueError(f"Unsupported Conda artifact: {artifact.name}")
    with zipfile.ZipFile(artifact) as archive:
        names = [
            name for name in archive.namelist()
            if (name.startswith("info-") or name.startswith("pkg-"))
            and name.endswith(".tar.zst")
        ]
        if len(names) != 2:
            raise ValueError(f"Expected info and payload archives in {artifact.name}")
        for name in names:
            with tempfile.NamedTemporaryFile(suffix=".tar") as uncompressed:
                subprocess.run(
                    ["zstd", "-dc"], input=archive.read(name),
                    stdout=uncompressed, check=True,
                )
                uncompressed.flush()
                with tarfile.open(uncompressed.name, "r:") as component:
                    extract_tar(component, destination)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")

    with args.inventory.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"name", "version", "build", "license", "binary_url", "sha256"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("Inventory has no rows or required columns")

    artifacts = args.output / "binary-artifacts"
    evidence = args.output / "package-evidence"
    manifest: list[dict[str, str]] = []
    for row in rows:
        filename = row["binary_url"].rsplit("/", 1)[-1]
        artifact = artifacts / filename
        download(row["binary_url"], artifact)
        if sha256(artifact) != row["sha256"].lower():
            raise SystemExit(f"Binary checksum mismatch: {row['binary_url']}")
        key = f"{row['name']}-{row['version']}-{row['build']}"
        package_evidence = evidence / key
        extract_evidence(artifact, package_evidence)
        files = sorted(
            path.relative_to(package_evidence).as_posix()
            for path in package_evidence.rglob("*")
            if path.is_file() and LICENSE_NAME.match(path.name)
        )
        about_path = package_evidence / "info" / "about.json"
        about = json.loads(about_path.read_text()) if about_path.is_file() else {}
        declared_files = about.get("license_file", "")
        if isinstance(declared_files, list):
            declared_files = " | ".join(str(item) for item in declared_files)
        manifest.append({
            "name": row["name"], "version": row["version"],
            "build": row["build"], "declared_license": row["license"],
            "binary_url": row["binary_url"], "binary_sha256": row["sha256"],
            "embedded_license_files": " | ".join(files),
            "recipe_declared_license_files": str(declared_files),
            "notice_status": "embedded" if files else "SOURCE_REVIEW_REQUIRED",
        })

    fields = list(manifest[0])
    with (args.output / "conda-notice-manifest.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    with (args.output / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in sorted(p for p in args.output.rglob("*") if p.is_file()
                           and p.name != "SHA256SUMS"):
            handle.write(f"{sha256(path)}  {path.relative_to(args.output)}\n")
    missing = sum(row["notice_status"] != "embedded" for row in manifest)
    print(
        f"Collected notice evidence for {len(manifest)} Conda packages; "
        f"{missing} require source-level notice review"
    )


if __name__ == "__main__":
    main()
