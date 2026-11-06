#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
WRAPPER="$PROJECT_DIR/scripts/run-cap-slurm.sh"
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

assembly="$tmp_dir/input genome.fasta"
outdir="$tmp_dir/results dir"
trash2="$tmp_dir/trash output"
mkdir -p "$trash2"
printf '>chr1\nACGT\n' > "$assembly"
printf 'repeat\n' > "$trash2/$(basename "$assembly")_repeats_with_seq.csv"
printf 'array\n' > "$trash2/$(basename "$assembly")_arrays.csv"

fake_conda="$tmp_dir/fake conda"
cat > "$fake_conda" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$@" > "${FAKE_CONDA_ARGS:?}"
exit "${FAKE_CONDA_STATUS:-0}"
EOF
chmod +x "$fake_conda"

expect_failure() {
    expected=$1
    shift
    if "$@" > "$tmp_dir/failure.out" 2>&1; then
        echo "Command unexpectedly succeeded: $*" >&2
        exit 1
    fi
    grep -F -- "$expected" "$tmp_dir/failure.out" >/dev/null || {
        echo "Expected error not found: $expected" >&2
        cat "$tmp_dir/failure.out" >&2
        exit 1
    }
}

expect_failure "submit this script with sbatch" \
    bash "$WRAPPER" --assembly "$assembly" --outdir "$outdir" \
    --conda "$fake_conda"

expect_failure "--outdir is required" \
    env CAP_SLURM_LOCAL_TEST=1 bash "$WRAPPER" \
    --assembly "$assembly" --conda "$fake_conda"

expect_failure "invalid --trash2 directory" \
    env CAP_SLURM_LOCAL_TEST=1 bash "$WRAPPER" \
    --assembly "$assembly" --outdir "$outdir" \
    --trash2 "$tmp_dir" --conda "$fake_conda"

expect_failure "--min-rep-size must not exceed --max-rep-size" \
    env CAP_SLURM_LOCAL_TEST=1 bash "$WRAPPER" \
    --assembly "$assembly" --outdir "$outdir" \
    --min-rep-size 20 --max-rep-size 10 --conda "$fake_conda"

args_file="$tmp_dir/args.txt"
run_tmp="$tmp_dir/node tmp"
mkdir -p "$run_tmp"
FAKE_CONDA_ARGS="$args_file" \
CAP_SLURM_LOCAL_TEST=1 \
SLURM_CPUS_PER_TASK=4 \
SLURM_JOB_ID=1234 \
SLURM_TMPDIR="$run_tmp" \
bash "$WRAPPER" \
    --assembly "$assembly" \
    --outdir "$outdir" \
    --trash2 "$trash2" \
    --max-rep-size 500 \
    --min-rep-size 8 \
    --env-name test-cap \
    --conda "$fake_conda" \
    > "$tmp_dir/success.out" 2>&1

python3 - "$args_file" "$assembly" "$outdir" "$trash2" "$run_tmp" <<'PY'
from pathlib import Path
import sys

args = Path(sys.argv[1]).read_text().splitlines()
work_dir = Path(args[args.index("-work-dir") + 1])
assert work_dir.parent == Path(sys.argv[5]).resolve(), work_dir
assert work_dir.name.startswith("cap-nextflow-1234-"), work_dir
assert len(work_dir.name.removeprefix("cap-nextflow-1234-")) == 8, work_dir
expected = [
    "run", "-n", "test-cap", "nextflow", "run", str(Path.cwd()),
    "--assembly", str(Path(sys.argv[2]).resolve()),
    "--cores", "4",
    "--outdir", str(Path(sys.argv[3]).resolve()),
    "-work-dir", str(work_dir),
    "--trash2", str(Path(sys.argv[4]).resolve()),
    "--max_rep_size", "500",
    "--min_rep_size", "8",
]
assert args == expected, (args, expected)
PY
work_dir=$(awk 'previous == "-work-dir" { print; exit } { previous=$0 }' "$args_file")
[ -d "$work_dir" ] || {
    echo "Successful run did not retain its work directory." >&2
    exit 1
}
grep -F -- "Nextflow work directory retained at: $work_dir" \
    "$tmp_dir/success.out" >/dev/null

without_trash_args="$tmp_dir/without-trash.txt"
FAKE_CONDA_ARGS="$without_trash_args" \
CAP_SLURM_LOCAL_TEST=1 \
SLURM_JOB_ID=1235 \
SLURM_TMPDIR="$run_tmp" \
bash "$WRAPPER" \
    --assembly "$assembly" \
    --outdir "$outdir" \
    --conda "$fake_conda" \
    > "$tmp_dir/without-trash.out" 2>&1
if grep -Fx -- "--trash2" "$without_trash_args" >/dev/null; then
    echo "Wrapper passed --trash2 when it was omitted." >&2
    exit 1
fi
second_work=$(awk 'previous == "-work-dir" { print; exit } { previous=$0 }' \
    "$without_trash_args")
[ "$second_work" != "$work_dir" ] || {
    echo "Separate runs reused the same work directory." >&2
    exit 1
}
[ -d "$second_work" ] || {
    echo "Second successful run did not retain its work directory." >&2
    exit 1
}

failed_output="$tmp_dir/failed.out"
if env FAKE_CONDA_ARGS="$tmp_dir/failed-args.txt" FAKE_CONDA_STATUS=9 \
    CAP_SLURM_LOCAL_TEST=1 SLURM_JOB_ID=1236 SLURM_TMPDIR="$run_tmp" \
    bash "$WRAPPER" --assembly "$assembly" --outdir "$outdir" \
    --conda "$fake_conda" > "$failed_output" 2>&1; then
    echo "Failed launcher command unexpectedly succeeded." >&2
    exit 1
fi
failed_work=$(awk -F ': ' '/Nextflow work directory retained at:/ { print $2 }' \
    "$failed_output")
case "$failed_work" in
    "$run_tmp"/cap-nextflow-1236-????????) ;;
    *) echo "Unexpected failed work directory: $failed_work" >&2; exit 1 ;;
esac
[ -d "$failed_work" ] || {
    echo "Failed run did not retain its work directory." >&2
    exit 1
}
grep -F -- "CAP failed; retain the work directory for diagnosis." \
    "$failed_output" >/dev/null

echo "Single-node SLURM wrapper tests passed without submitting a job."
