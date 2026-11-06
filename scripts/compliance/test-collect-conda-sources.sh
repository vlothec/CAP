#!/usr/bin/env bash

set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
mkdir -p "$tmpdir/http" "$tmpdir/package/info/recipe"

printf 'exact upstream source\n' > "$tmpdir/http/source-1.0.tar.gz"
source_md5="$(md5sum "$tmpdir/http/source-1.0.tar.gz" | cut -d' ' -f1)"
cat > "$tmpdir/package/info/recipe/meta.yaml" <<EOF
package:
  name: example
  version: 1.0
source:
  - url:
      - http://127.0.0.1:8765/missing-source.tar.gz
      - http://127.0.0.1:8765/source-1.0.tar.gz
    md5: $source_md5
build:
  number: 0
EOF
cat > "$tmpdir/package/info/about.json" <<'EOF'
{"extra":{"remote_url":"https://example.invalid/feedstock","sha":"abc123"}}
EOF
cat > "$tmpdir/package/info/index.json" <<'EOF'
{"name":"example","version":"1.0","build":"0"}
EOF
tar -cjf "$tmpdir/http/example-1.0-0.tar.bz2" -C "$tmpdir/package" info
binary_sha="$(sha256sum "$tmpdir/http/example-1.0-0.tar.bz2" | cut -d' ' -f1)"
cat > "$tmpdir/review.csv" <<EOF
name,version,build,license,review_class,collect_source,review_note,\
binary_url,binary_sha256
example,1.0,0,GPL-3.0-only,gpl-source-review,yes,test,\
http://127.0.0.1:8765/example-1.0-0.tar.bz2,$binary_sha
EOF

python3 -m http.server 8765 --directory "$tmpdir/http" \
    >"$tmpdir/server.log" 2>&1 &
server_pid=$!
trap 'kill "$server_pid" 2>/dev/null || true; rm -rf "$tmpdir"' EXIT
sleep 1
"$PROJECT_DIR/scripts/compliance/collect-conda-sources.py" \
    "$tmpdir/review.csv" "$tmpdir/output"
grep -F 'https://example.invalid/feedstock,abc123' \
    "$tmpdir/output/source-to-binary.csv" >/dev/null
grep -F "upstream-source/example-1.0-0/01-source-1.0.tar.gz" \
    "$tmpdir/output/source-to-binary.csv" >/dev/null
(cd "$tmpdir/output" && sha256sum --check --strict SHA256SUMS >/dev/null)

echo "Conda source collector test passed."
