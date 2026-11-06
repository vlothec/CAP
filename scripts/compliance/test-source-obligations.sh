#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

cat > "$tmpdir/inventory.csv" <<'CSV'
name,version,build,license,binary_url,sha256
gpl,1,0,GPL-3.0-only,https://example.invalid/gpl,aaa
lgpl,1,0,LGPL-2.1-or-later,https://example.invalid/lgpl,bbb
alternative,1,0,BSD-2-Clause OR GPL-2.0-or-later,https://example.invalid/alt,ccc
exception,1,0,GPL-2.0 WITH Classpath-exception-2.0,https://example.invalid/exc,ddd
permissive,1,0,MIT,https://example.invalid/mit,eee
CSV

"$PROJECT_DIR/scripts/compliance/classify-source-obligations.py" \
    "$tmpdir/inventory.csv" "$tmpdir/output.csv"

grep -F 'gpl-source-review,REVIEW' "$tmpdir/output.csv" >/dev/null
grep -F 'lgpl-source-review,REVIEW' "$tmpdir/output.csv" >/dev/null
grep -F 'alternative-license-review,REVIEW' "$tmpdir/output.csv" >/dev/null
grep -F 'exception-review,REVIEW' "$tmpdir/output.csv" >/dev/null
! grep -F 'permissive' "$tmpdir/output.csv" >/dev/null
test "$(wc -l < "$tmpdir/output.csv")" -eq 5

echo "Source-obligation classifier test passed."
