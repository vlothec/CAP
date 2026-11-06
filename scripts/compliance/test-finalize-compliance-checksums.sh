#!/usr/bin/env bash

set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

mkdir -p "$tmpdir/bundle/nested"
printf 'root\n' > "$tmpdir/bundle/root.txt"
printf 'nested\n' > "$tmpdir/bundle/nested/data.txt"
printf 'obsolete\n' > "$tmpdir/bundle/nested/SHA256SUMS"
printf 'obsolete\n' > "$tmpdir/bundle/nested/data.txt.sha256"

"$PROJECT_DIR/scripts/compliance/finalize-compliance-checksums.sh" \
    "$tmpdir/bundle" >/dev/null

test -f "$tmpdir/bundle/SHA256SUMS"
test ! -e "$tmpdir/bundle/nested/SHA256SUMS"
test ! -e "$tmpdir/bundle/nested/data.txt.sha256"
test "$(find "$tmpdir/bundle" -type f -name SHA256SUMS | wc -l)" -eq 1
test "$(find "$tmpdir/bundle" -type f -name '*.sha256' | wc -l)" -eq 0
(
    cd "$tmpdir/bundle"
    sha256sum --check --strict SHA256SUMS >/dev/null
)

echo "Compliance checksum finalizer test passed."
