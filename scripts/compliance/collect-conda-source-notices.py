#!/usr/bin/env python3

"""Extract source-level license evidence from verified Conda source archives."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
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


def wanted(name: str) -> bool:
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts:
        return False
    return LICENSE_NAME.match(path.name) is not None or path.name == "DESCRIPTION"


def safe_target(root: Path, archive_name: str) -> Path:
    path = PurePosixPath(archive_name)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Unsafe archive path: {archive_name}")
    return root.joinpath(*path.parts)


def extract_archive(archive: Path, destination: Path) -> list[str]:
    files: list[str] = []
    try:
        with tarfile.open(archive, "r:*") as source:
            for member in source:
                if not member.isfile() or not wanted(member.name):
                    continue
                target = safe_target(destination, member.name)
                target.parent.mkdir(parents=True, exist_ok=True)
                extracted = source.extractfile(member)
                if extracted is None:
                    continue
                with target.open("wb") as output:
                    shutil.copyfileobj(extracted, output)
                files.append(target.relative_to(destination).as_posix())
        return files
    except tarfile.ReadError:
        pass
    try:
        with zipfile.ZipFile(archive) as source:
            for name in source.namelist():
                if name.endswith("/") or not wanted(name):
                    continue
                target = safe_target(destination, name)
                target.parent.mkdir(parents=True, exist_ok=True)
                with source.open(name) as input_file, target.open("wb") as output:
                    shutil.copyfileobj(input_file, output)
                files.append(target.relative_to(destination).as_posix())
        return files
    except zipfile.BadZipFile:
        return []


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("notice_manifest", type=Path)
    parser.add_argument("source_mapping", type=Path)
    parser.add_argument("source_bundle", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--shared-license-dir", action="append", type=Path, default=[])
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")

    with args.notice_manifest.open(newline="", encoding="utf-8") as handle:
        notices = list(csv.DictReader(handle))
    with args.source_mapping.open(newline="", encoding="utf-8") as handle:
        mappings = list(csv.DictReader(handle))
    gaps = {
        (row["name"], row["version"], row["build"]): row
        for row in notices if row["notice_status"] != "embedded"
    }
    by_package: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in mappings:
        key = (row["package"], row["version"], row["build"])
        by_package.setdefault(key, []).append(row)

    args.output.mkdir(parents=True)
    shared = args.output / "shared-licenses"
    for directory in args.shared_license_dir:
        shutil.copytree(directory, shared / directory.name)

    records: list[dict[str, str]] = []
    for key, notice in sorted(gaps.items()):
        package_dir = args.output / "package-evidence" / "-".join(key)
        files: list[str] = []
        sources = by_package.get(key, [])
        for index, source in enumerate(sources, start=1):
            archive = args.source_bundle / source["collected_path"]
            archive_dir = package_dir / f"source-{index:02d}"
            extracted = extract_archive(archive, archive_dir)
            files.extend(
                f"source-{index:02d}/{name}" for name in extracted
            )
        license_files = [
            name for name in files if LICENSE_NAME.match(PurePosixPath(name).name)
        ]
        description_files = [
            name for name in files if PurePosixPath(name).name == "DESCRIPTION"
        ]
        status = "source-license-files" if license_files else "REVIEW_REQUIRED"
        records.append({
            "name": key[0], "version": key[1], "build": key[2],
            "declared_license": notice["declared_license"],
            "source_archives": " | ".join(
                row["collected_path"] for row in sources
            ),
            "source_license_files": " | ".join(license_files),
            "source_description_files": " | ".join(description_files),
            "review_status": status,
        })

    fields = list(records[0])
    with (args.output / "source-notice-manifest.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    with (args.output / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in sorted(p for p in args.output.rglob("*") if p.is_file()
                           and p.name != "SHA256SUMS"):
            handle.write(f"{sha256(path)}  {path.relative_to(args.output)}\n")
    unresolved = sum(row["review_status"] == "REVIEW_REQUIRED" for row in records)
    print(
        f"Extracted source notice evidence for {len(records)} packages; "
        f"{unresolved} require manual review"
    )


if __name__ == "__main__":
    main()
