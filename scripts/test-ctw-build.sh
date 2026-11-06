#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

fake_conda="$TMP_DIR/conda"
cat > "$fake_conda" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
shift 3
exec "$@"
EOF
chmod +x "$fake_conda"

fixture="$TMP_DIR/project"
mkdir -p "$fixture/bin/src/BCT"
cat > "$fixture/bin/src/BCT/Makefile" <<'EOF'
all:
	cp ctw-calc.sh ../../ctw-calc
	chmod +x ../../ctw-calc
EOF

cat > "$fixture/bin/src/BCT/ctw-calc.sh" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
cat > "$2" <<'CSV'
chromosome,bin_mid,bin_value
smoke,61.00,20.57
CSV
EOF

CAP_PROJECT_DIR="$fixture" CAP_CONDA_EXE="$fake_conda" \
    CAP_ENV_NAME=unused bash "$PROJECT_DIR/scripts/build-ctw.sh" \
    > "$TMP_DIR/success.log"
grep -Fq "ctw-calc compiled and passed its execution smoke test." \
    "$TMP_DIR/success.log"

cat > "$fixture/bin/src/BCT/ctw-calc.sh" <<'EOF'
#!/usr/bin/env bash
exit 126
EOF

set +e
CAP_PROJECT_DIR="$fixture" CAP_CONDA_EXE="$fake_conda" \
    CAP_ENV_NAME=unused bash "$PROJECT_DIR/scripts/build-ctw.sh" \
    > "$TMP_DIR/failure.log" 2>&1
status=$?
set -e

if [ "$status" -eq 0 ]; then
    echo "CTW build validation accepted a binary that could not execute." >&2
    exit 1
fi

echo "CTW build validation tests passed."
