#!/usr/bin/env python3

"""Collect exact Conda recipe evidence and checksum-verified upstream source."""

from __future__ import annotations

import argparse
import bz2
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
import zipfile

HASH_RE = re.compile(r"^[0-9a-fA-F]{64}$")
SOURCE_HASH_RE = re.compile(
    r"^\s*(sha256|sha512|md5):\s*([0-9a-fA-F]{32,128})\s*$"
)
URL_RE = re.compile(r"^\s*(?:-\s*)?url:\s*(.*?)\s*$")
TARGET_RE = re.compile(r"^\s*target_directory:\s*(\S.*?)\s*$")
SOURCE_RE = re.compile(r"^(\s*)source:\s*$")
TOP_KEY_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_-]*:\s*")


def file_hash(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256(path: Path) -> str:
    return file_hash(path)


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "CAP-compliance/1"})
    with urllib.request.urlopen(request) as response, destination.open("wb") as output:
        shutil.copyfileobj(response, output)


def safe_member(name: str) -> bool:
    path = Path(name)
    return not path.is_absolute() and ".." not in path.parts


def extract_info(artifact: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    if artifact.name.endswith(".tar.bz2"):
        with tarfile.open(artifact, "r:bz2") as archive:
            members = [m for m in archive if safe_member(m.name) and (
                m.name == "info/about.json" or m.name == "info/index.json" or
                m.name.startswith("info/recipe/") or
                m.name.startswith("info/licenses/")
            )]
            archive.extractall(destination, members=members, filter="data")
        return
    if not artifact.name.endswith(".conda"):
        raise ValueError(f"Unsupported Conda artifact: {artifact.name}")
    with zipfile.ZipFile(artifact) as archive:
        info_names = [name for name in archive.namelist()
                      if name.startswith("info-") and name.endswith(".tar.zst")]
        if len(info_names) != 1:
            raise ValueError(f"Expected one info archive in {artifact.name}")
        with tempfile.NamedTemporaryFile(suffix=".tar") as uncompressed:
            process = subprocess.run(
                ["zstd", "-dc"], input=archive.read(info_names[0]),
                stdout=uncompressed, check=True,
            )
            del process
            uncompressed.flush()
            with tarfile.open(uncompressed.name, "r:") as info_archive:
                members = [m for m in info_archive if safe_member(m.name) and (
                    m.name == "info/about.json" or m.name == "info/index.json" or
                    m.name.startswith("info/recipe/") or
                    m.name.startswith("info/licenses/")
                )]
                info_archive.extractall(destination, members=members, filter="data")


def rendered_recipe(recipe_dir: Path) -> Path:
    rendered = recipe_dir / "rendered_recipe.yaml"
    if rendered.is_file():
        return rendered
    legacy = recipe_dir / "meta.yaml"
    if legacy.is_file() and "{{" not in legacy.read_text(errors="replace"):
        return legacy
    raise ValueError(f"No fully rendered recipe in {recipe_dir}")


def source_block(lines: list[str]) -> list[str]:
    candidates: list[tuple[int, int]] = []
    for index, line in enumerate(lines):
        match = SOURCE_RE.match(line)
        if match:
            candidates.append((index, len(match.group(1))))
    for start, indent in candidates:
        block: list[str] = []
        for line in lines[start + 1:]:
            line_indent = len(line) - len(line.lstrip())
            if (line.strip() and line_indent <= indent and
                    not line.lstrip().startswith("-")):
                break
            block.append(line)
        if any(URL_RE.match(line) is not None for line in block):
            return block
    raise ValueError("Rendered recipe has no source block containing a URL")


def clean_scalar(value: str) -> str:
    value = value.strip().strip("'\"")
    if "{{" in value or "${{" in value:
        raise ValueError(f"Unresolved source template: {value}")
    return value


def parse_sources(recipe: Path) -> list[dict[str, str]]:
    block = source_block(recipe.read_text(errors="strict").splitlines())
    sources: list[dict[str, str]] = []
    current: dict[str, str] = {}
    urls: list[str] = []
    reading_urls = False

    def finish() -> None:
        nonlocal current, urls, reading_urls
        if not urls and current.get("url"):
            urls = [current["url"]]
        if urls or current:
            if not urls or "hash_value" not in current:
                raise ValueError(f"Incomplete source mapping in {recipe}")
            current["url"] = urls[0]
            if len(urls) > 1:
                current["urls"] = "\n".join(urls)
            sources.append(current)
        current = {}
        urls = []
        reading_urls = False

    for line in block:
        stripped = line.strip()
        indent = len(line) - len(line.lstrip())
        url = URL_RE.match(line)
        digest = SOURCE_HASH_RE.match(line)
        target = TARGET_RE.match(line)
        if stripped.startswith("- url:") and (urls or current):
            finish()
            url = URL_RE.match(line)
        if reading_urls and stripped.startswith("- http"):
            urls.append(clean_scalar(stripped[2:]))
            continue
        if reading_urls and stripped and not stripped.startswith("#"):
            reading_urls = False
        if url:
            value = clean_scalar(url.group(1))
            if value:
                urls.append(value)
            else:
                reading_urls = True
        elif digest:
            current["hash_algorithm"] = digest.group(1).lower()
            current["hash_value"] = digest.group(2).lower()
        elif target:
            current["target_directory"] = clean_scalar(target.group(1))
    finish()
    if not sources:
        raise ValueError(f"No source mappings parsed from {recipe}")
    return sources


def artifact_name(url: str) -> str:
    name = Path(urllib.request.url2pathname(url.split("?", 1)[0])).name
    if not name:
        raise ValueError(f"Source URL has no filename: {url}")
    return name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("review_csv", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--conservative", action="store_true",
                        help="collect every candidate, including REVIEW rows")
    args = parser.parse_args()

    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")
    with args.review_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"name", "version", "build", "collect_source", "binary_url",
                "binary_sha256"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit("Review CSV has no rows or required columns")
    selected = [row for row in rows if row["collect_source"].lower() == "yes" or
                (args.conservative and row["collect_source"] == "REVIEW")]
    if not selected:
        raise SystemExit("No rows selected for source collection")

    args.output.mkdir(parents=True)
    artifacts = args.output / "binary-artifacts"
    evidence = args.output / "package-evidence"
    sources_dir = args.output / "upstream-source"
    mappings: list[dict[str, str]] = []

    for row in selected:
        binary_url = row["binary_url"]
        binary = artifacts / artifact_name(binary_url)
        download(binary_url, binary)
        expected_binary = row["binary_sha256"].lower()
        if not HASH_RE.fullmatch(expected_binary) or sha256(binary) != expected_binary:
            raise SystemExit(f"Binary checksum mismatch: {binary_url}")
        package_key = f"{row['name']}-{row['version']}-{row['build']}"
        package_evidence = evidence / package_key
        extract_info(binary, package_evidence)
        recipe = rendered_recipe(package_evidence / "info" / "recipe")
        about_path = package_evidence / "info" / "about.json"
        about = json.loads(about_path.read_text()) if about_path.is_file() else {}
        recipe_sha = about.get("extra", {}).get("sha", "")
        recipe_url = about.get("extra", {}).get("remote_url", "")
        for number, source in enumerate(parse_sources(recipe), start=1):
            source_urls = source.get("urls", source["url"]).splitlines()
            source_url = source_urls[0]
            source_hash = source["hash_value"]
            source_hash_algorithm = source["hash_algorithm"]
            package_source_dir = sources_dir / package_key
            partial = package_source_dir / f"{number:02d}.partial"
            errors = []
            for candidate_url in source_urls:
                try:
                    download(candidate_url, partial)
                    source_url = candidate_url
                    source_name = artifact_name(source_url)
                    destination = package_source_dir / f"{number:02d}-{source_name}"
                    partial.replace(destination)
                    break
                except Exception as error:
                    partial.unlink(missing_ok=True)
                    errors.append(f"{candidate_url}: {error}")
            else:
                raise SystemExit("No source mirror succeeded:\n" + "\n".join(errors))
            actual_recipe_hash = file_hash(destination, source_hash_algorithm)
            if actual_recipe_hash != source_hash:
                raise SystemExit(f"Source checksum mismatch: {source_url}")
            mappings.append({
                "package": row["name"], "version": row["version"],
                "build": row["build"], "binary_url": binary_url,
                "binary_sha256": expected_binary, "recipe_repository": recipe_url,
                "recipe_commit": recipe_sha, "source_url": source_url,
                "recipe_hash_algorithm": source_hash_algorithm,
                "recipe_hash_value": source_hash,
                "source_sha256": sha256(destination),
                "target_directory": source.get("target_directory", ""),
                "collected_path": str(destination.relative_to(args.output)),
                "review_status": row["collect_source"],
            })

    fields = list(mappings[0])
    with (args.output / "source-to-binary.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(mappings)
    with (args.output / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in sorted(p for p in args.output.rglob("*") if p.is_file() and
                           p.name != "SHA256SUMS"):
            handle.write(f"{sha256(path)}  {path.relative_to(args.output)}\n")
    print(f"Collected {len(mappings)} source mappings for {len(selected)} packages")


if __name__ == "__main__":
    main()
