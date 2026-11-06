#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR=${CAP_PROJECT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}
ENV_NAME=${CAP_ENV_NAME:-cap-pipeline}
CONDA_EXE=${CAP_CONDA_EXE:-conda}

"$CONDA_EXE" run -n "$ENV_NAME" make -C "$PROJECT_DIR/bin/src/BCT"

work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT

python3 - "$work_dir/input.fa" <<'PY'
import sys

with open(sys.argv[1], "w", encoding="ascii") as fasta:
    fasta.write(">smoke\n" + "ACGT" * 30 + "\n")
PY

"$PROJECT_DIR/bin/ctw-calc" \
    "$work_dir/input.fa" "$work_dir/output.csv" 100 >/dev/null

cat > "$work_dir/expected.csv" <<'EOF'
chromosome,bin_mid,bin_value
smoke,61.00,20.57
EOF

if ! cmp -s "$work_dir/expected.csv" "$work_dir/output.csv"; then
    echo "ctw-calc smoke-test output did not match the expected result." >&2
    diff -u "$work_dir/expected.csv" "$work_dir/output.csv" >&2 || true
    exit 1
fi

echo "ctw-calc compiled and passed its execution smoke test."
