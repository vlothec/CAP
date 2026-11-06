#!/usr/bin/env bash

set -euo pipefail

IMAGE="${1:-cap:local}"
OUTDIR="${2:-oci-inventory}"
RUNTIME="${CONTAINER_RUNTIME:-podman}"

command -v "$RUNTIME" >/dev/null || {
    echo "Container runtime not found: $RUNTIME" >&2
    exit 1
}
mkdir -p "$OUTDIR"

image_id="$($RUNTIME image inspect --format '{{.Id}}' "$IMAGE")"
image_digest="$($RUNTIME image inspect --format '{{.Digest}}' "$IMAGE")"
image_architecture="$($RUNTIME image inspect --format '{{.Architecture}}' "$IMAGE")"
image_os="$($RUNTIME image inspect --format '{{.Os}}' "$IMAGE")"
python3 - "$image_id" "$image_digest" "$image_architecture" "$image_os" \
    > "$OUTDIR/image-metadata.json" <<'PY'
import json
import sys

keys = ("Id", "Digest", "Architecture", "Os")
json.dump(dict(zip(keys, sys.argv[1:])), sys.stdout, indent=2)
print()
PY

"$RUNTIME" run --rm --entrypoint /bin/bash "$IMAGE" -c '
set -euo pipefail
printf "name,version,build,channel,license\n"
python3 - <<"PY"
import csv
import glob
import json
import sys

writer = csv.writer(sys.stdout, lineterminator="\n")
for path in sorted(glob.glob("/opt/cap-env/conda-meta/*.json")):
    with open(path, encoding="utf-8") as handle:
        item = json.load(handle)
    writer.writerow([
        item.get("name", ""),
        item.get("version", ""),
        item.get("build", ""),
        item.get("channel", ""),
        item.get("license", "UNKNOWN") or "UNKNOWN",
    ])
PY
' > "$OUTDIR/conda-packages.csv"

"$RUNTIME" run --rm --entrypoint /bin/bash "$IMAGE" -c '
set -euo pipefail
printf "package,version\n"
dpkg-query -W -f="\${Package},\${Version}\\n" | LC_ALL=C sort
' > "$OUTDIR/debian-packages.csv"

"$RUNTIME" run --rm --entrypoint /bin/bash "$IMAGE" -c '
set -euo pipefail
printf "component,revision_or_hash,license_expression,license_or_notice_file\n"
printf "CAP,%s,%s,%s\n" \
    "$(sha256sum /opt/CAP/main.nf | cut -d" " -f1)" \
    "MIT" \
    "/opt/CAP/LICENSE"
printf "BCT_ctw-calc,%s,%s,%s\n" \
    "$(sha256sum /opt/CAP/bin/ctw-calc | cut -d" " -f1)" \
    "GPL-2.0-or-later" \
    "/opt/CAP/bin/src/BCT/README.md"
printf "TRASH_2,%s,%s,%s\n" \
    "$(sha256sum /opt/CAP/modules/TRASH_2/src/TRASH.R | cut -d" " -f1)" \
    "MIT" \
    "/opt/CAP/modules/TRASH_2/license.txt"
printf "model,%s,%s,%s\n" \
    "$(sha256sum /opt/CAP/model/centromeric_model_v2.pkl | cut -d" " -f1)" \
    "NOASSERTION" \
    "/opt/CAP/model/centromeric_model_v2.pkl"
' > "$OUTDIR/embedded-components.csv"

(
    cd "$OUTDIR"
    sha256sum image-metadata.json conda-packages.csv debian-packages.csv \
        embedded-components.csv > SHA256SUMS
)

echo "OCI inventory written to $OUTDIR"
