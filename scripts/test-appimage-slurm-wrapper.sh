#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
WRAPPER="$PROJECT_DIR/scripts/run-cap-appimage-slurm.sh"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

assembly="$TMP_DIR/input genome.fasta"
outdir="$TMP_DIR/results dir"
trash2="$TMP_DIR/trash output"
te_gff="$TMP_DIR/TE annotation.gff3"
gene_gff="$TMP_DIR/gene annotation.gff"
metadata="$TMP_DIR/sample metadata.csv"
templates="$TMP_DIR/repeat templates.fasta"
mkdir -p "$trash2"
printf '>chr1\nACGT\n' > "$assembly"
printf 'repeat\n' > "$trash2/$(basename "$assembly")_repeats_with_seq.csv"
printf 'array\n' > "$trash2/$(basename "$assembly")_arrays.csv"
printf '##gff-version 3\n' > "$te_gff"
printf '##gff-version 3\n' > "$gene_gff"
printf 'chromosome,length\nchr1,4\n' > "$metadata"
printf '>repeat\nACGT\n' > "$templates"

fake_appimage="$TMP_DIR/CAP test.AppImage"
cat > "$fake_appimage" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$@" > "${FAKE_APPIMAGE_ARGS:?}"
exit "${FAKE_APPIMAGE_STATUS:-0}"
EOF
chmod +x "$fake_appimage"

expect_failure() {
    local expected=$1
    shift
    if "$@" > "$TMP_DIR/failure.out" 2>&1; then
        echo "Command unexpectedly succeeded: $*" >&2
        exit 1
    fi
    grep -F -- "$expected" "$TMP_DIR/failure.out" >/dev/null || {
        echo "Expected error not found: $expected" >&2
        cat "$TMP_DIR/failure.out" >&2
        exit 1
    }
}

expect_failure "submit this script with sbatch" \
    bash "$WRAPPER" --appimage "$fake_appimage" \
    --assembly "$assembly" --outdir "$outdir"
expect_failure "--appimage is required" \
    env CAP_SLURM_LOCAL_TEST=1 bash "$WRAPPER" \
    --assembly "$assembly" --outdir "$outdir"
expect_failure "invalid --trash2 directory" \
    env CAP_SLURM_LOCAL_TEST=1 bash "$WRAPPER" \
    --appimage "$fake_appimage" --assembly "$assembly" --outdir "$outdir" \
    --trash2 "$TMP_DIR"
expect_failure "--min-rep-size must not exceed --max-rep-size" \
    env CAP_SLURM_LOCAL_TEST=1 bash "$WRAPPER" \
    --appimage "$fake_appimage" --assembly "$assembly" --outdir "$outdir" \
    --min-rep-size 20 --max-rep-size 10
[ ! -e "$outdir" ] || {
    echo "Invalid input created the output directory." >&2
    exit 1
}
expect_failure "SLURM_CPUS_PER_TASK must be a positive integer" \
    env CAP_SLURM_LOCAL_TEST=1 SLURM_CPUS_PER_TASK=zero bash "$WRAPPER" \
    --appimage "$fake_appimage" --assembly "$assembly" --outdir "$outdir"
[ ! -e "$outdir" ] || {
    echo "Invalid CPU input created the output directory." >&2
    exit 1
}

args_file="$TMP_DIR/args.txt"
node_tmp="$TMP_DIR/node tmp"
mkdir -p "$node_tmp"
FAKE_APPIMAGE_ARGS="$args_file" \
CAP_SLURM_LOCAL_TEST=1 \
SLURM_CPUS_PER_TASK=6 \
SLURM_JOB_ID=4321 \
SLURM_TMPDIR="$node_tmp" \
bash "$WRAPPER" \
    --appimage "$fake_appimage" \
    --assembly "$assembly" \
    --outdir "$outdir" \
    --trash2 "$trash2" \
    --te-gff "$te_gff" \
    --gene-gff "$gene_gff" \
    --metadata "$metadata" \
    --templates "$templates" \
    --max-rep-size 500 \
    --min-rep-size 8 > "$TMP_DIR/success.out" 2>&1

python3 - "$args_file" "$assembly" "$outdir" "$trash2" "$te_gff" \
    "$gene_gff" "$metadata" "$templates" "$node_tmp" <<'PY'
from pathlib import Path
import sys

args = Path(sys.argv[1]).read_text().splitlines()
work_dir = Path(args[args.index("-work-dir") + 1])
assert work_dir.parent == Path(sys.argv[9]).resolve(), work_dir
assert work_dir.name.startswith("cap-appimage-nextflow-4321-"), work_dir
assert len(work_dir.name.removeprefix("cap-appimage-nextflow-4321-")) == 8, work_dir
expected = [
    "--assembly", str(Path(sys.argv[2]).resolve()),
    "--cores", "6",
    "--outdir", str(Path(sys.argv[3]).resolve()),
    "-work-dir", str(work_dir),
    "--trash2", str(Path(sys.argv[4]).resolve()),
    "--te_gff", str(Path(sys.argv[5]).resolve()),
    "--gene_gff", str(Path(sys.argv[6]).resolve()),
    "--metadata", str(Path(sys.argv[7]).resolve()),
    "--templates", str(Path(sys.argv[8]).resolve()),
    "--max_rep_size", "500",
    "--min_rep_size", "8",
]
assert args == expected, (args, expected)
PY
success_work=$(awk 'previous == "-work-dir" { print; exit } { previous=$0 }' \
    "$args_file")
[ -d "$success_work" ] || {
    echo "Successful run did not retain its work directory." >&2
    exit 1
}
grep -F -- "Nextflow work directory retained at: $success_work" \
    "$TMP_DIR/success.out" >/dev/null

without_optional="$TMP_DIR/without-optional.txt"
FAKE_APPIMAGE_ARGS="$without_optional" \
CAP_SLURM_LOCAL_TEST=1 \
SLURM_JOB_ID=4322 \
SLURM_TMPDIR="$node_tmp" \
bash "$WRAPPER" --appimage "$fake_appimage" \
    --assembly "$assembly" --outdir "$outdir" > "$TMP_DIR/minimal.out" 2>&1
if grep -Eq -- '--trash2|--te_gff|--gene_gff|--metadata|--templates' \
    "$without_optional"; then
    echo "Wrapper passed an omitted optional argument." >&2
    exit 1
fi
second_work=$(awk 'previous == "-work-dir" { print; exit } { previous=$0 }' \
    "$without_optional")
[ "$second_work" != "$success_work" ] || {
    echo "Separate runs reused the same work directory." >&2
    exit 1
}
[ -d "$second_work" ] || {
    echo "Second successful run did not retain its work directory." >&2
    exit 1
}

failed_output="$TMP_DIR/failed.out"
if env FAKE_APPIMAGE_ARGS="$TMP_DIR/failed-args.txt" FAKE_APPIMAGE_STATUS=9 \
    CAP_SLURM_LOCAL_TEST=1 SLURM_JOB_ID=4323 SLURM_TMPDIR="$node_tmp" \
    bash "$WRAPPER" --appimage "$fake_appimage" \
    --assembly "$assembly" --outdir "$outdir" > "$failed_output" 2>&1; then
    echo "Failed AppImage launcher command unexpectedly succeeded." >&2
    exit 1
fi
failed_work=$(awk -F ': ' '/Nextflow work directory retained at:/ { print $2 }' \
    "$failed_output")
case "$failed_work" in
    "$node_tmp"/cap-appimage-nextflow-4323-????????) ;;
    *) echo "Unexpected failed work directory: $failed_work" >&2; exit 1 ;;
esac
[ -d "$failed_work" ] || {
    echo "Failed run did not retain its work directory." >&2
    exit 1
}
grep -F -- "CAP AppImage failed; retain the work directory for diagnosis." \
    "$failed_output" >/dev/null

printf 'AppImage SLURM wrapper tests passed without submitting a job.\n'
