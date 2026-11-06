#!/usr/bin/env bash

set -euo pipefail
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
mkdir -p "$tmpdir/http" "$tmpdir/package/info/recipe" \
    "$tmpdir/package/share/example"
printf 'recipe license\n' > "$tmpdir/package/info/recipe/LICENSE.txt"
printf 'payload notice\n' > "$tmpdir/package/share/example/NOTICE.md"
printf '%s\n' '{"license_file":"LICENSE.txt"}' \
    > "$tmpdir/package/info/about.json"
printf '%s\n' '{"name":"example"}' > "$tmpdir/package/info/index.json"
tar -cjf "$tmpdir/http/example-1.0-0.tar.bz2" -C "$tmpdir/package" .
binary_sha="$(sha256sum "$tmpdir/http/example-1.0-0.tar.bz2" | cut -d' ' -f1)"
cat > "$tmpdir/inventory.csv" <<EOF
name,version,build,license,binary_url,sha256
example,1.0,0,MIT,http://127.0.0.1:8766/example-1.0-0.tar.bz2,$binary_sha
EOF
python3 -m http.server 8766 --directory "$tmpdir/http" \
    >"$tmpdir/server.log" 2>&1 &
server_pid=$!
trap 'kill "$server_pid" 2>/dev/null || true; rm -rf "$tmpdir"' EXIT
sleep 1
"$PROJECT_DIR/scripts/compliance/collect-conda-notices.py" \
    "$tmpdir/inventory.csv" "$tmpdir/output"
grep -F 'info/recipe/LICENSE.txt | share/example/NOTICE.md' \
    "$tmpdir/output/conda-notice-manifest.csv" >/dev/null
(cd "$tmpdir/output" && sha256sum --check --strict SHA256SUMS >/dev/null)
echo "Conda notice collector test passed."
