#!/usr/bin/env bash

set -euo pipefail

BUNDLE="${1:?Usage: finalize-compliance-checksums.sh BUNDLE}"
test -d "$BUNDLE" || {
    echo "Bundle directory not found: $BUNDLE" >&2
    exit 1
}

find "$BUNDLE" -mindepth 2 -type f \
    \( -name SHA256SUMS -o -name '*.sha256' \) -delete
rm -f "$BUNDLE/SHA256SUMS"

(
    cd "$BUNDLE"
    find . -type f ! -name SHA256SUMS -print0 | LC_ALL=C sort -z | \
        xargs -0 sha256sum > SHA256SUMS
)

echo "Finalized bundle checksum manifest: $BUNDLE/SHA256SUMS"
