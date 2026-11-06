#!/usr/bin/env bash

set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
mkdir -p "$tmpdir/inventory"
cat > "$tmpdir/inventory/conda-packages.csv" <<EOF
name,version,build,channel,license
example,1.0,h123_0,https://conda.example/linux-64,MIT
example-r,2.0,r44_1,https://conda.example/noarch,GPL-3.0-only
EOF
cat > "$tmpdir/inventory/debian-packages.csv" <<EOF
package,version
base-files,12.4
libc6,2.36-9
EOF
a_hash="$(printf 'a%.0s' {1..64})"
b_hash="$(printf 'b%.0s' {1..64})"
c_hash="$(printf 'c%.0s' {1..64})"
d_hash="$(printf 'd%.0s' {1..64})"
cat > "$tmpdir/inventory/embedded-components.csv" <<EOF
component,revision_or_hash,license_expression,license_or_notice_file
CAP,$a_hash,MIT,/opt/CAP/LICENSE
BCT_ctw-calc,$b_hash,GPL-2.0-or-later,/opt/CAP/bin/src/BCT/README.md
TRASH_2,$c_hash,MIT,/opt/CAP/modules/TRASH_2/license.txt
model,$d_hash,NOASSERTION,/opt/CAP/model/centromeric_model_v2.pkl
EOF
cat > "$tmpdir/inventory/image-metadata.json" <<EOF
{
  "Id": "0123", "Digest": "sha256:4567",
  "Architecture": "amd64", "Os": "linux"
}
EOF
"$PROJECT_DIR/scripts/compliance/generate-cyclonedx-sbom.py" \
    "$tmpdir/inventory" "$tmpdir/sbom.json"
python3 - "$tmpdir/sbom.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    document = json.load(handle)

assert document["bomFormat"] == "CycloneDX"
assert document["specVersion"] == "1.6"
assert document["metadata"]["component"]["name"] == "CAP"
components = {component["name"]: component for component in document["components"]}
assert set(components) == {
    "BCT_ctw-calc", "TRASH_2", "base-files", "example", "example-r",
    "libc6", "model",
}
assert components["example"]["purl"] == (
    "pkg:conda/example@1.0?build=h123_0&channel="
    "https%3A%2F%2Fconda.example%2Flinux-64"
)
assert components["base-files"]["purl"] == "pkg:deb/debian/base-files@12.4"
assert components["BCT_ctw-calc"]["licenses"] == [
    {"expression": "GPL-2.0-or-later"}
]
assert components["example-r"]["licenses"] == [{"expression": "GPL-3.0-only"}]
assert components["model"]["licenses"] == [
    {"license": {"name": "NOASSERTION"}}
]
root_ref = document["metadata"]["component"]["bom-ref"]
assert document["dependencies"] == [{
    "ref": root_ref,
    "dependsOn": [component["bom-ref"] for component in document["components"]],
}]
PY
echo "CycloneDX SBOM generation test passed."
