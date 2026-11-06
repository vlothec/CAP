#!/usr/bin/env python3

"""Generate a deterministic CycloneDX SBOM from CAP's exact OCI inventory."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import urllib.parse
import uuid

NAMESPACE = uuid.UUID("ce891ff3-e49c-50eb-bd88-6dbaab3f76af")
SPDX_EXPRESSION = re.compile(
    r"^[A-Za-z0-9.+-]+(?: WITH [A-Za-z0-9.+-]+)?"
    r"(?: (?:AND|OR) [A-Za-z0-9.+-]+(?: WITH [A-Za-z0-9.+-]+)?)*$"
)
SPDX_IDS = {
    "0BSD", "Apache-2.0", "Artistic-1.0-Perl", "Artistic-2.0",
    "BSD-2-Clause", "BSD-3-Clause", "Classpath-exception-2.0", "FTL",
    "GCC-exception-3.1", "GPL-1.0-or-later", "GPL-2.0-only",
    "GPL-2.0-or-later", "GPL-3.0-only", "GPL-3.0-or-later", "HPND",
    "IJG", "ISC", "LGPL-2.0-or-later", "LGPL-2.1-only",
    "LGPL-2.1-or-later", "LGPL-3.0-only", "LGPL-3.0-or-later", "MIT",
    "MPL-1.1", "OFL-1.1", "Python-2.0", "TCL", "X11", "Zlib",
    "bzip2-1.0.6", "curl", "zlib-acknowledgement",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def encoded(value: str) -> str:
    return urllib.parse.quote(value, safe=".-_~")


def bom_ref(kind: str, *values: str) -> str:
    return f"urn:uuid:{uuid.uuid5(NAMESPACE, '|'.join((kind, *values)))}"


def valid_spdx_expression(value: str) -> bool:
    if not SPDX_EXPRESSION.fullmatch(value):
        return False
    identifiers = re.split(r" (?:AND|OR|WITH) ", value)
    return all(item in SPDX_IDS or item.startswith("LicenseRef-")
               for item in identifiers)


def license_choice(value: str) -> dict[str, dict[str, str]]:
    value = value.strip() or "UNKNOWN"
    if valid_spdx_expression(value) and "provisional" not in value.lower():
        return {"expression": value}
    return {"license": {"name": value}}


def sha256_hash(value: str) -> list[dict[str, str]]:
    if re.fullmatch(r"[0-9a-fA-F]{64}", value):
        return [{"alg": "SHA-256", "content": value.lower()}]
    return []


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("inventory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")

    conda = read_csv(args.inventory / "conda-packages.csv")
    debian = read_csv(args.inventory / "debian-packages.csv")
    embedded = read_csv(args.inventory / "embedded-components.csv")
    image_data = json.loads(
        (args.inventory / "image-metadata.json").read_text(encoding="utf-8")
    )
    if len(conda) != len({(r["name"], r["version"], r["build"]) for r in conda}):
        raise SystemExit("Duplicate exact Conda package record")
    if len(debian) != len({(r["package"], r["version"]) for r in debian}):
        raise SystemExit("Duplicate exact Debian package record")

    root_row = next((row for row in embedded if row["component"] == "CAP"), None)
    if root_row is None:
        raise SystemExit("Embedded component inventory has no CAP row")
    root_ref = bom_ref("embedded", "CAP", root_row["revision_or_hash"])
    root = {
        "type": "application",
        "bom-ref": root_ref,
        "name": "CAP",
        "version": root_row["revision_or_hash"][:12],
        "hashes": sha256_hash(root_row["revision_or_hash"]),
        "licenses": [license_choice(root_row["license_expression"])],
        "properties": [
            {"name": "cap:notice-file", "value": root_row["license_or_notice_file"]},
            {"name": "oci:image:id", "value": image_data.get("Id", "")},
            {"name": "oci:image:digest", "value": image_data.get("Digest", "")},
            {
                "name": "oci:image:architecture",
                "value": image_data.get("Architecture", ""),
            },
            {"name": "oci:image:os", "value": image_data.get("Os", "")},
        ],
    }

    components: list[dict[str, object]] = []
    for row in sorted(conda, key=lambda r: (r["name"], r["version"], r["build"])):
        ref = bom_ref("conda", row["name"], row["version"], row["build"])
        qualifiers = urllib.parse.urlencode({
            "build": row["build"], "channel": row["channel"],
        })
        components.append({
            "type": "library", "bom-ref": ref, "name": row["name"],
            "version": row["version"],
            "purl": (
                f"pkg:conda/{encoded(row['name'])}@{encoded(row['version'])}"
                f"?{qualifiers}"
            ),
            "licenses": [license_choice(row["license"])],
            "properties": [
                {"name": "conda:build", "value": row["build"]},
                {"name": "conda:channel", "value": row["channel"]},
            ],
        })
    for row in sorted(debian, key=lambda r: (r["package"], r["version"])):
        ref = bom_ref("debian", row["package"], row["version"])
        components.append({
            "type": "library", "bom-ref": ref, "name": row["package"],
            "version": row["version"],
            "purl": (
                f"pkg:deb/debian/{encoded(row['package'])}"
                f"@{encoded(row['version'])}"
            ),
            "properties": [{"name": "cap:package-manager", "value": "dpkg"}],
        })
    for row in sorted(embedded, key=lambda r: r["component"]):
        if row["component"] == "CAP":
            continue
        ref = bom_ref("embedded", row["component"], row["revision_or_hash"])
        component_type = "data" if row["component"] == "model" else "application"
        components.append({
            "type": component_type, "bom-ref": ref,
            "name": row["component"],
            "version": row["revision_or_hash"][:12],
            "hashes": sha256_hash(row["revision_or_hash"]),
            "licenses": [license_choice(row["license_expression"])],
            "properties": [{
                "name": "cap:notice-file",
                "value": row["license_or_notice_file"],
            }],
        })

    image_key = image_data.get("Digest") or image_data.get("Id") or root_ref
    timestamp = image_data.get("Created")
    try:
        datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
    except (AttributeError, ValueError):
        timestamp = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    document = {
        "bomFormat": "CycloneDX", "specVersion": "1.6", "version": 1,
        "serialNumber": f"urn:uuid:{uuid.uuid5(NAMESPACE, image_key)}",
        "metadata": {
            "timestamp": timestamp,
            "tools": {"components": [{
                "type": "application", "name": "CAP inventory SBOM generator",
                "version": "1",
            }]},
            "component": root,
        },
        "components": components,
        "dependencies": [{
            "ref": root_ref,
            "dependsOn": [component["bom-ref"] for component in components],
        }],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    checksum = hashlib.sha256(args.output.read_bytes()).hexdigest()
    print(
        f"Generated CycloneDX SBOM with {len(conda)} Conda, {len(debian)} Debian, "
        f"and {len(embedded)} embedded components; SHA-256 {checksum}"
    )


if __name__ == "__main__":
    main()
