#!/usr/bin/env bash

set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
mkdir -p "$tmpdir/source/upstream-source/example-1.0-0/pkg" \
    "$tmpdir/shared"
printf 'Package: example\nVersion: 1.0\nLicense: GPL-2\n' \
    > "$tmpdir/source/upstream-source/example-1.0-0/pkg/DESCRIPTION"
tar -czf "$tmpdir/source/upstream-source/example-1.0-0/source.tar.gz" \
    -C "$tmpdir/source/upstream-source/example-1.0-0" pkg
rm -rf "$tmpdir/source/upstream-source/example-1.0-0/pkg"
cat > "$tmpdir/notices.csv" <<EOF
name,version,build,declared_license,notice_status
example,1.0,0,GPL-2,SOURCE_REVIEW_REQUIRED
metadata-only,2.0,0,MIT,SOURCE_REVIEW_REQUIRED
EOF
cat > "$tmpdir/mapping.csv" <<EOF
package,version,build,collected_path
example,1.0,0,upstream-source/example-1.0-0/source.tar.gz
EOF
printf 'canonical GPL text\n' > "$tmpdir/shared/GPL-2"
"$PROJECT_DIR/scripts/compliance/collect-conda-source-notices.py" \
    "$tmpdir/notices.csv" "$tmpdir/mapping.csv" "$tmpdir/source" \
    "$tmpdir/evidence" --shared-license-dir "$tmpdir/shared"
cat > "$tmpdir/policy.csv" <<EOF
name,version,build,resolution,canonical_license_files,review_note
example,1.0,0,R-shared-license,GPL-2,DESCRIPTION refers to canonical GPL-2.
metadata-only,2.0,0,recipe-license-declaration,,Exact recipe is retained as evidence.
EOF
"$PROJECT_DIR/scripts/compliance/finalize-conda-source-notices.py" \
    "$tmpdir/evidence/source-notice-manifest.csv" "$tmpdir/policy.csv" \
    "$tmpdir/shared" "$tmpdir/final.csv"
grep -F 'R-shared-license,GPL-2' "$tmpdir/final.csv" >/dev/null
grep -F 'recipe-license-declaration' "$tmpdir/final.csv" >/dev/null
(cd "$tmpdir/evidence" && sha256sum --check --strict SHA256SUMS >/dev/null)
echo "Conda source-notice review test passed."
