#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
LAUNCHER="$PROJECT_DIR/scripts/run-trash2.sh"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

BIN_DIR="$TMP_DIR/bin"
LOG_DIR="$TMP_DIR/logs"
mkdir -p "$BIN_DIR" "$LOG_DIR"
printf '>sequence\nACGTACGTACGTACGT\n' > "$TMP_DIR/input genome.fasta"
printf '>template\nACGTACGT\n' > "$TMP_DIR/template repeats.fasta"

cat > "$BIN_DIR/conda" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
: "${FAKE_CONDA_LOG:?}"
printf '%q ' "$@" >> "$FAKE_CONDA_LOG"
printf '\n' >> "$FAKE_CONDA_LOG"
if [ "${1:-}" = env ] && [ "${2:-}" = list ]; then
    printf '# conda environments:\n#\n                         /tmp/test-cap\n'
    exit 0
fi
[ "${1:-}" = run ] || exit 90
shift
[ "${1:-}" = --no-capture-output ] || exit 91
shift
[ "${1:-}" = -p ] || exit 92
shift
[ "${1:-}" = /tmp/test-cap ] || exit 93
shift
exec "$@"
EOF

cat > "$BIN_DIR/Rscript" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
: "${FAKE_RSCRIPT_LOG:?}"
printf 'R_REPOSITORIES=%s\n' "${R_REPOSITORIES:-}" >> "$FAKE_RSCRIPT_LOG"
printf 'R_PROFILE_USER=%s\n' "${R_PROFILE_USER:-}" >> "$FAKE_RSCRIPT_LOG"
[ -s "${R_PROFILE_USER:-}" ] || exit 97
if [ -n "${CAPTURED_R_PROFILE:-}" ] && [ ! -e "$CAPTURED_R_PROFILE" ]; then
    cp "$R_PROFILE_USER" "$CAPTURED_R_PROFILE"
fi
[ "${http_proxy:-}" = 'http://127.0.0.1:9' ] || exit 99
[ "${https_proxy:-}" = 'http://127.0.0.1:9' ] || exit 100
printf '%q ' "$@" >> "$FAKE_RSCRIPT_LOG"
printf '\n' >> "$FAKE_RSCRIPT_LOG"
if [ "${1:-}" = -e ]; then
    [ "${FAKE_PREFLIGHT_FAIL:-0}" != 1 ] || exit 42
    exit 0
fi
script=$1
shift
[ "$script" = "${EXPECTED_TRASH2_SCRIPT:?}" ] || exit 95
assembly=
outdir=
templates=
cores=
max_rep=
min_rep=
while [ "$#" -gt 0 ]; do
    case "$1" in
        -f) assembly=$2; shift 2 ;;
        -o) outdir=$2; shift 2 ;;
        --cores_no) cores=$2; shift 2 ;;
        --max_rep_size) max_rep=$2; shift 2 ;;
        --min_rep_size) min_rep=$2; shift 2 ;;
        --templates) templates=$2; shift 2 ;;
        *) exit 96 ;;
    esac
done
[ -f "$assembly" ]
[ -d "$outdir" ]
[ -f "$templates" ]
[ "$cores" = 2 ]
[ "$max_rep" = 300 ]
[ "$min_rep" = 8 ]
printf 'repeat output\n' > "$outdir/$(basename "$assembly")_repeats_with_seq.csv"
printf 'array output\n' > "$outdir/$(basename "$assembly")_arrays.csv"
EOF
chmod +x "$BIN_DIR/conda" "$BIN_DIR/Rscript"

run_launcher() {
    PATH="$BIN_DIR:$PATH" \
    FAKE_CONDA_LOG="$LOG_DIR/conda.log" \
    FAKE_RSCRIPT_LOG="$LOG_DIR/rscript.log" \
    CAPTURED_R_PROFILE="$TMP_DIR/generated.Rprofile" \
    EXPECTED_TRASH2_SCRIPT="$PROJECT_DIR/modules/TRASH_2/src/TRASH.R" \
        "$LAUNCHER" "$@"
}

expect_failure() {
    local description=$1
    shift
    if run_launcher "$@" >"$TMP_DIR/failure.log" 2>&1; then
        echo "Expected failure: $description" >&2
        exit 1
    fi
}

run_launcher \
    --assembly "$TMP_DIR/input genome.fasta" \
    --outdir "$TMP_DIR/TRASH_2 output" \
    --env-name test-cap \
    --conda "$BIN_DIR/conda" \
    --cores 2 \
    --max-rep-size 300 \
    --min-rep-size 8 \
    --templates "$TMP_DIR/template repeats.fasta"

test -f "$TMP_DIR/TRASH_2 output/input genome.fasta_repeats_with_seq.csv"
test -f "$TMP_DIR/TRASH_2 output/input genome.fasta_arrays.csv"

REAL_RSCRIPT=${CAP_TEST_RSCRIPT:-}
R_RUNNER=()
if [ -n "$REAL_RSCRIPT" ]; then
    R_RUNNER=("$REAL_RSCRIPT")
elif command -v Rscript >/dev/null 2>&1; then
    R_RUNNER=(Rscript)
elif command -v "${CAP_TEST_CONDA:-conda}" >/dev/null 2>&1; then
    R_RUNNER=("${CAP_TEST_CONDA:-conda}" run --no-capture-output
        -n "${CAP_TEST_ENV_NAME:-cap-pipeline}" Rscript)
fi
if [ "${#R_RUNNER[@]}" -gt 0 ]; then
    if R_PROFILE_USER="$TMP_DIR/generated.Rprofile" \
        "${R_RUNNER[@]}" -e 'install.packages("cap-runtime-install-must-fail")' \
        >"$TMP_DIR/install-guard.log" 2>&1; then
        echo "Runtime R package installation was not blocked." >&2
        exit 1
    fi
    grep -F 'runtime package installation is disabled' \
        "$TMP_DIR/install-guard.log" >/dev/null || {
        echo "R profile failed for an unexpected reason." >&2
        cat "$TMP_DIR/install-guard.log" >&2
        exit 1
    }
    echo "Runtime R package installation guard passed."
else
    echo "Skipping R installation-guard behavior test: no Rscript or Conda available."
fi

test "$(grep -c '^R_REPOSITORIES=' "$LOG_DIR/rscript.log")" -eq 2
test "$(grep -c '^R_PROFILE_USER=' "$LOG_DIR/rscript.log")" -eq 2
if grep '^R_REPOSITORIES=' "$LOG_DIR/rscript.log" |
    grep -vFx 'R_REPOSITORIES=@CRAN@=file:///nonexistent/cap-no-runtime-installs' >/dev/null; then
    echo "Launcher did not force offline R repositories." >&2
    exit 1
fi
python3 - "$LOG_DIR/conda.log" <<'PY'
from pathlib import Path
import shlex
import sys
runs = [shlex.split(line) for line in Path(sys.argv[1]).read_text().splitlines()]
conda_runs = [run for run in runs if run and run[0] == "run"]
assert len(conda_runs) == 2, conda_runs
for run in conda_runs:
    assert run[:4] == ["run", "--no-capture-output", "-p", "/tmp/test-cap"], run
PY

expect_failure "missing assembly" --outdir "$TMP_DIR/missing-assembly"
expect_failure "missing output directory" --assembly "$TMP_DIR/input genome.fasta"
expect_failure "invalid core count" \
    --assembly "$TMP_DIR/input genome.fasta" --outdir "$TMP_DIR/bad-cores" --cores 0
expect_failure "minimum larger than maximum" \
    --assembly "$TMP_DIR/input genome.fasta" --outdir "$TMP_DIR/bad-range" \
    --min-rep-size 20 --max-rep-size 10
expect_failure "missing template" \
    --assembly "$TMP_DIR/input genome.fasta" --outdir "$TMP_DIR/bad-template" \
    --templates "$TMP_DIR/not-there.fasta"
expect_failure "unknown option" \
    --assembly "$TMP_DIR/input genome.fasta" --outdir "$TMP_DIR/unknown" --bogus value

rm -f "$LOG_DIR/conda.log" "$LOG_DIR/rscript.log"
if PATH="$BIN_DIR:$PATH" \
    FAKE_CONDA_LOG="$LOG_DIR/conda.log" \
    FAKE_RSCRIPT_LOG="$LOG_DIR/rscript.log" \
    FAKE_PREFLIGHT_FAIL=1 \
    EXPECTED_TRASH2_SCRIPT="$PROJECT_DIR/modules/TRASH_2/src/TRASH.R" \
    "$LAUNCHER" \
        --assembly "$TMP_DIR/input genome.fasta" \
        --outdir "$TMP_DIR/preflight-failure" \
        --env-name test-cap >"$TMP_DIR/preflight.log" 2>&1; then
    echo "Expected dependency preflight failure." >&2
    exit 1
fi
test ! -e "$TMP_DIR/preflight-failure"
test "$(grep -c '^R_REPOSITORIES=' "$LOG_DIR/rscript.log")" -eq 1
test "$(grep -c '^R_PROFILE_USER=' "$LOG_DIR/rscript.log")" -eq 1

printf 'Standalone TRASH_2 launcher tests passed.\n'
