#!/usr/bin/env bash

set -euo pipefail

IMAGE="${1:-cap:local}"
OUTDIR="${2:-oci-compliance-bundle}"
RUNTIME="${CONTAINER_RUNTIME:-podman}"
GIT_COMMAND="${GIT_COMMAND:-git}"
INCLUDE_OCI_ARCHIVE="${INCLUDE_OCI_ARCHIVE:-0}"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

command -v "$RUNTIME" >/dev/null || {
    echo "Container runtime not found: $RUNTIME" >&2
    exit 1
}
git_available=0
if [[ "$GIT_COMMAND" == */* ]]; then
    [[ -x "$GIT_COMMAND" ]] && git_available=1
elif command -v "$GIT_COMMAND" >/dev/null; then
    git_available=1
fi
if [[ "$git_available" != "1" ]]; then
    echo "Git is required to validate the compliance source checkout." >&2
    exit 1
fi

test ! -e "$OUTDIR" || {
    echo "Output path already exists: $OUTDIR" >&2
    exit 1
}
mkdir -p "$OUTDIR/inventory" "$OUTDIR/licenses" "$OUTDIR/source" \
    "$OUTDIR/sbom"

"$PROJECT_DIR/scripts/compliance/generate-oci-inventory.sh" \
    "$IMAGE" "$OUTDIR/inventory"

image_id="$($RUNTIME image inspect --format '{{.Id}}' "$IMAGE")"
image_digest="$($RUNTIME image inspect --format '{{index .RepoDigests 0}}' \
    "$IMAGE" 2>/dev/null || true)"
image_cap_revision="$($RUNTIME image inspect --format \
    '{{index .Config.Labels "org.opencontainers.image.revision"}}' \
    "$IMAGE" 2>/dev/null || true)"
image_trash_revision="$($RUNTIME image inspect --format \
    '{{index .Config.Labels "org.cap-project.trash-2.revision"}}' \
    "$IMAGE" 2>/dev/null || true)"
cap_revision="$image_cap_revision"
trash_revision="$image_trash_revision"
if [[ "$git_available" == "1" ]]; then
    if [[ -n "$($GIT_COMMAND -C "$PROJECT_DIR" status --porcelain \
        --untracked-files=normal --ignore-submodules=none)" ]]; then
        echo "Compliance bundles require a clean checkout with no untracked files." >&2
        exit 1
    fi
    checkout_cap_revision="$($GIT_COMMAND -C "$PROJECT_DIR" rev-parse HEAD)"
    checkout_trash_revision="$($GIT_COMMAND -C "$PROJECT_DIR" \
        rev-parse HEAD:modules/TRASH_2)"
    mapfile -t resolved_revisions < <(
        "$PROJECT_DIR/scripts/compliance/resolve-release-revisions.sh" \
            "$image_cap_revision" "$image_trash_revision" \
            "$checkout_cap_revision" "$checkout_trash_revision"
    )
    [[ "${#resolved_revisions[@]}" == "2" ]] || {
        echo "Could not resolve CAP and TRASH_2 release revisions." >&2
        exit 1
    }
    cap_revision="${resolved_revisions[0]}"
    trash_revision="${resolved_revisions[1]}"
fi
{
    printf 'image_reference=%s\n' "$IMAGE"
    printf 'image_id=%s\n' "$image_id"
    printf 'registry_digest=%s\n' "${image_digest:-UNAVAILABLE_FOR_LOCAL_IMAGE}"
    printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'cap_revision=%s\n' "${cap_revision:-UNAVAILABLE}"
    printf 'trash_2_revision=%s\n' "${trash_revision:-UNAVAILABLE}"
} > "$OUTDIR/RELEASE-METADATA.txt"

cp "$PROJECT_DIR/LICENSE" "$OUTDIR/licenses/CAP-MIT.txt"
cp "$PROJECT_DIR/THIRD_PARTY_NOTICES.md" "$OUTDIR/licenses/"
cp "$PROJECT_DIR/docs/DEPENDENCY_SOURCE_REVIEW.md" "$OUTDIR/"
cp "$PROJECT_DIR/licenses/GPL-2.0.txt" "$OUTDIR/licenses/"
cp "$PROJECT_DIR/licenses/GPL-3.0.txt" "$OUTDIR/licenses/"
cp "$PROJECT_DIR/modules/TRASH_2/license.txt" \
    "$OUTDIR/licenses/TRASH_2-MIT.txt"

tar --sort=name --mtime='UTC 1970-01-01' --owner=0 --group=0 \
    --numeric-owner -C "$PROJECT_DIR" -cJf \
    "$OUTDIR/source/BCT-corresponding-source.tar.xz" bin/src/BCT

"$RUNTIME" run --rm --interactive --entrypoint python3 "$IMAGE" - <<'PY' \
    > "$OUTDIR/inventory/conda-source-urls.csv"
import csv
import glob
import json
import sys

writer = csv.writer(sys.stdout, lineterminator="\n")
writer.writerow(["name", "version", "build", "license", "binary_url", "sha256"])
for path in sorted(glob.glob("/opt/cap-env/conda-meta/*.json")):
    with open(path, encoding="utf-8") as handle:
        item = json.load(handle)
    writer.writerow([
        item.get("name", ""),
        item.get("version", ""),
        item.get("build", ""),
        item.get("license", "UNKNOWN") or "UNKNOWN",
        item.get("url", ""),
        item.get("sha256", ""),
    ])
PY

test -s "$OUTDIR/inventory/conda-source-urls.csv" || {
    echo "Conda source URL inventory is empty." >&2
    exit 1
}
grep -q '^name,version,build,license,binary_url,sha256' \
    "$OUTDIR/inventory/conda-source-urls.csv"
"$PROJECT_DIR/scripts/compliance/classify-source-obligations.py" \
    "$OUTDIR/inventory/conda-source-urls.csv" \
    "$OUTDIR/source-obligations.csv"

"$PROJECT_DIR/scripts/compliance/generate-cyclonedx-sbom.py" \
    "$OUTDIR/inventory" "$OUTDIR/sbom/cap-primary.cdx.json"

if [[ "$INCLUDE_OCI_ARCHIVE" == "1" ]]; then
    runtime_name="${RUNTIME##*/}"
    if [[ "$runtime_name" != "podman" ]]; then
        echo "INCLUDE_OCI_ARCHIVE=1 requires Podman; Docker save creates a" \
            "Docker archive, not an OCI archive." >&2
        exit 1
    fi
    command -v zstd >/dev/null || {
        echo "zstd is required when INCLUDE_OCI_ARCHIVE=1." >&2
        exit 1
    }
    archive="$OUTDIR/image.oci.tar.zst"
    "$RUNTIME" save --format oci-archive "$IMAGE" | \
        zstd -T0 -19 -o "$archive"
fi

cat > "$OUTDIR/README.txt" <<'EOF'
CAP OCI compliance bundle
=========================

This bundle records the image identity, package inventories, declared licenses,
repository revisions, and complete corresponding source for CAP's BCT-derived
ctw-calc executable. Package binary URLs identify exact Conda artifacts and are
a starting point for corresponding-source collection; they are not themselves
source archives. Apply
scripts/compliance/apply-source-review-policy.py to source-obligations.csv,
then run scripts/compliance/collect-conda-sources.py against the reviewed CSV.
Collect embedded package notices with
scripts/compliance/collect-conda-notices.py. For packages flagged for source
review, collect the reviewed supplemental inputs in
licenses/conda-notice-supplemental-sources.csv with
scripts/compliance/collect-conda-notice-sources.py, extract source evidence with
scripts/compliance/collect-conda-source-notices.py, and apply the exact reviewed
policy in licenses/conda-source-notice-review.csv with
scripts/compliance/finalize-conda-source-notices.py. Export and collect Debian
copyright/source evidence with
scripts/compliance/export-debian-compliance-inputs.sh and
scripts/compliance/collect-debian-sources.py. These network-dependent
collections are external bundle steps, not part of basic bundle generation. The
--conservative option preserves all Conda candidates while review is pending;
it does not replace legal review. Obtain and preserve exact corresponding
source for packages whose licenses require it before public distribution.

An OCI archive is included only when INCLUDE_OCI_ARCHIVE=1 and the container
runtime is Podman. Docker users can create a separate Docker archive with
`docker image save` or convert the daemon image to OCI format with Skopeo. The
primary CycloneDX SBOM is generated directly from the exact image inventories.
EOF

"$PROJECT_DIR/scripts/compliance/finalize-compliance-checksums.sh" "$OUTDIR"

echo "OCI compliance bundle written to $OUTDIR"
