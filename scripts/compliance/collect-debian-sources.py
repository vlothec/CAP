#!/usr/bin/env python3

"""Collect exact Debian source archives for an OCI runtime package inventory."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import lzma
from pathlib import Path
import shutil
import tarfile
import urllib.parse
import urllib.request

INDEX_URLS = (
    "https://deb.debian.org/debian/dists/bookworm/main/source/Sources.xz",
    "https://deb.debian.org/debian/dists/bookworm-updates/main/source/Sources.xz",
    "https://deb.debian.org/debian-security/dists/bookworm-security/main/source/"
    "Sources.xz",
)
SNAPSHOT = "https://snapshot.debian.org"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "CAP-compliance/1"})
    with urllib.request.urlopen(request) as response:
        return response.read()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(fetch(url))


def parse_stanzas(text: str) -> list[dict[str, str]]:
    stanzas: list[dict[str, str]] = []
    for block in text.split("\n\n"):
        item: dict[str, str] = {}
        key = ""
        for line in block.splitlines():
            if line.startswith(" ") and key:
                item[key] += "\n" + line[1:]
            elif ":" in line and not line.startswith(" "):
                key, value = line.split(":", 1)
                item[key] = value.lstrip()
        if item:
            stanzas.append(item)
    return stanzas


def current_sources() -> dict[tuple[str, str], dict[str, str]]:
    result: dict[tuple[str, str], dict[str, str]] = {}
    for url in INDEX_URLS:
        text = lzma.decompress(fetch(url)).decode()
        for item in parse_stanzas(text):
            key = (item.get("Package", ""), item.get("Version", ""))
            if all(key):
                result[key] = item
    return result


def stanza_files(item: dict[str, str]) -> list[dict[str, str]]:
    files = []
    for line in item.get("Checksums-Sha256", "").splitlines():
        if not line.strip():
            continue
        digest, size, name = line.split(maxsplit=2)
        files.append({"sha256": digest, "size": size, "name": name})
    if not files:
        raise ValueError(
            "No SHA-256 source files for "
            f"{item.get('Package')} {item.get('Version')}"
        )
    return files


def snapshot_files(package: str, version: str) -> list[dict[str, str]]:
    quoted_package = urllib.parse.quote(package, safe="")
    quoted_version = urllib.parse.quote(version, safe="")
    url = f"{SNAPSHOT}/mr/package/{quoted_package}/{quoted_version}/srcfiles"
    hashes = json.loads(fetch(url))["result"]
    files = []
    for entry in hashes:
        digest = entry["hash"]
        info_url = f"{SNAPSHOT}/mr/file/{digest}/info"
        records = json.loads(fetch(info_url))["result"]
        preferred = next(
            (
                record for record in records
                if record["archive_name"] == "debian-security"
            ),
            records[0],
        )
        files.append({
            "sha1": digest, "name": preferred["name"],
            "size": str(preferred["size"]),
            "url": f"{SNAPSHOT}/file/{digest}",
        })
    return sorted(files, key=lambda item: item["name"])


def safe_extract(archive_path: Path, destination: Path) -> None:
    with tarfile.open(archive_path, "r:gz") as archive:
        archive.extractall(destination, filter="data")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("binary_source_map", type=Path)
    parser.add_argument("copyright_archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")

    with args.binary_source_map.open(newline="", encoding="utf-8") as handle:
        mappings = list(csv.DictReader(handle))
    required = {
        "binary_package", "binary_version", "source_package", "source_version",
    }
    if not mappings or not required.issubset(mappings[0]):
        raise SystemExit("Binary/source map has no rows or required columns")

    args.output.mkdir(parents=True)
    notices = args.output / "installed-copyright"
    notices.mkdir()
    safe_extract(args.copyright_archive, notices)
    shutil.copy2(args.binary_source_map, args.output / "binary-to-source.csv")

    index = current_sources()
    source_pairs = sorted({
        (row["source_package"], row["source_version"]) for row in mappings
    })
    records: list[dict[str, str]] = []
    for package, version in source_pairs:
        destination = args.output / "source" / f"{package}-{version}"
        item = index.get((package, version))
        if item:
            directory = item["Directory"]
            for source_file in stanza_files(item):
                url = (
                    f"https://deb.debian.org/debian/{directory}/"
                    f"{source_file['name']}"
                )
                if directory.startswith("pool/updates/"):
                    url = (
                        f"https://deb.debian.org/debian-security/{directory}/"
                        f"{source_file['name']}"
                    )
                path = destination / source_file["name"]
                download(url, path)
                if sha256(path) != source_file["sha256"]:
                    raise SystemExit(f"Debian source checksum mismatch: {url}")
                records.append({
                    "source_package": package, "source_version": version,
                    "filename": source_file["name"], "download_url": url,
                    "sha256": source_file["sha256"], "evidence": "Sources index",
                })
        else:
            for source_file in snapshot_files(package, version):
                path = destination / source_file["name"]
                download(source_file["url"], path)
                records.append({
                    "source_package": package, "source_version": version,
                    "filename": source_file["name"],
                    "download_url": source_file["url"], "sha256": sha256(path),
                    "evidence": f"Debian Snapshot SHA-1 {source_file['sha1']}",
                })

    fields = list(records[0])
    with (args.output / "source-files.csv").open(
            "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    with (args.output / "SHA256SUMS").open("w", encoding="utf-8") as handle:
        for path in sorted(p for p in args.output.rglob("*") if p.is_file()
                           and p.name != "SHA256SUMS"):
            handle.write(f"{sha256(path)}  {path.relative_to(args.output)}\n")
    print(
        f"Collected {len(records)} source files for {len(source_pairs)} "
        f"Debian source packages"
    )


if __name__ == "__main__":
    main()
