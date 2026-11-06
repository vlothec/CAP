#!/usr/bin/env bash

set -euo pipefail

ENV_NAME="${CAP_ENV_NAME:-cap-pipeline}"
CONDA_EXE="${CAP_CONDA_EXE:-conda}"
PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
ASSEMBLY="${ASSEMBLY:-test/ath_Chr1_extraction_trc.fasta}"
RESULTS_DIR="${RESULTS_DIR:-results_test_conda}"
WORK_DIR="${WORK_DIR:-work_test_conda}"

if [ ! -f "$ASSEMBLY" ]; then
    echo "Assembly FASTA not found: $ASSEMBLY" >&2
    exit 1
fi

if [ ! -f modules/TRASH_2/src/TRASH.R ]; then
    echo "TRASH_2 submodule is not initialized; run make install first." >&2
    exit 1
fi

if [ ! -x bin/ctw-calc ]; then
    echo "bin/ctw-calc is missing; run make install or make compile first." >&2
    exit 1
fi

if ! "$CONDA_EXE" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "Conda environment '$ENV_NAME' does not exist; run make install first." >&2
    exit 1
fi

rm -rf "$RESULTS_DIR" "$WORK_DIR"

"$CONDA_EXE" run -n "$ENV_NAME" nextflow run . \
    --assembly "$ASSEMBLY" \
    --cores 1 \
    --outdir "$RESULTS_DIR" \
    -work-dir "$WORK_DIR"

trash2_command=$(find "$WORK_DIR" -type f -name .command.sh -exec \
    grep -l -- 'modules/TRASH_2/src/TRASH.R' {} + | head -n 1)
if [ -z "$trash2_command" ] || \
    ! grep -F -- '--min_rep_size 7' "$trash2_command" >/dev/null; then
    echo "TRASH_2 task did not receive the configured minimum repeat size." >&2
    exit 1
fi

file_count=$(find "$RESULTS_DIR" -maxdepth 1 -type f | wc -l)
if [ "$file_count" -ne 16 ]; then
    echo "Expected 16 result files, found $file_count in $RESULTS_DIR." >&2
    exit 1
fi

for extension in csv png rds; do
    if ! find "$RESULTS_DIR" -maxdepth 1 -type f -name "*.$extension" -print -quit | grep -q .; then
        echo "No .$extension result found in $RESULTS_DIR." >&2
        exit 1
    fi
done

if find "$RESULTS_DIR" -maxdepth 1 -type f -empty -print -quit | grep -q .; then
    echo "An empty result file was produced in $RESULTS_DIR." >&2
    exit 1
fi

EXPECTED_CHECKSUMS="$PROJECT_DIR/test/expected-results.sha256"
if [ ! -f "$EXPECTED_CHECKSUMS" ]; then
    echo "Expected checksum manifest not found: $EXPECTED_CHECKSUMS" >&2
    exit 1
fi

if ! (cd "$RESULTS_DIR" && sha256sum --check "$EXPECTED_CHECKSUMS"); then
    echo "Result checksums do not match $EXPECTED_CHECKSUMS." >&2
    exit 1
fi

echo "Conda smoke test passed: 16 nonempty files match expected checksums."
