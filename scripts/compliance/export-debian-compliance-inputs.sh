#!/usr/bin/env bash

set -euo pipefail
IMAGE="${1:?Usage: export-debian-compliance-inputs.sh IMAGE OUTDIR}"
OUTDIR="${2:?Usage: export-debian-compliance-inputs.sh IMAGE OUTDIR}"
RUNTIME="${CONTAINER_RUNTIME:-podman}"

test ! -e "$OUTDIR" || {
    echo "Output path already exists: $OUTDIR" >&2
    exit 1
}
command -v "$RUNTIME" >/dev/null || {
    echo "Container runtime not found: $RUNTIME" >&2
    exit 1
}
mkdir -p "$OUTDIR"

"$RUNTIME" run --rm --entrypoint /bin/bash "$IMAGE" -c '
format="\${binary:Package},\${Version},\${source:Package},\${source:Version}\n"
dpkg-query -W -f="$format" |
    LC_ALL=C sort
' | {
    printf '%s\n' \
        'binary_package,binary_version,source_package,source_version'
    cat
} > "$OUTDIR/binary-to-source.csv"

"$RUNTIME" run --rm --entrypoint /bin/bash "$IMAGE" -c '
set -euo pipefail
cd /
find usr/share/doc -mindepth 2 -maxdepth 2 -type f -name copyright \
    -print0 | LC_ALL=C sort -z |
    tar --null --files-from=- --create --gzip --file=-
' > "$OUTDIR/installed-copyright.tar.gz"

(
    cd "$OUTDIR"
    sha256sum binary-to-source.csv installed-copyright.tar.gz > SHA256SUMS
)
echo "Debian compliance inputs written to $OUTDIR"
