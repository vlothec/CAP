#!/usr/bin/env bash

set -euo pipefail

CAP_DIR=/opt/CAP
RESULTS_DIR=/work/container-test-results
WORK_DIR=/work/container-test-work
EXPECTED="$CAP_DIR/test/expected-results.sha256"

for command in nextflow java Rscript python3 mafft nhmmer; do
    command -v "$command" >/dev/null || {
        echo "Required command missing from image: $command" >&2
        exit 1
    }
done

if command -v conda >/dev/null; then
    echo "Conda package manager must not be present in the runtime image." >&2
    exit 1
fi
for build_command in gcc g++ gfortran make cmake git; do
    if command -v "$build_command" >/dev/null; then
        echo "Build command must not be present in runtime image: $build_command" >&2
        exit 1
    fi
done

test -x "$CAP_DIR/bin/ctw-calc"
test -f "$CAP_DIR/modules/TRASH_2/src/TRASH.R"
test -f "$CAP_DIR/modules/TRASH_2/dep/HOR.V3.3"
test ! -e "$CAP_DIR/modules/TRASH_2/dep/hmmer"
test ! -e "$CAP_DIR/modules/TRASH_2/dep/mafft-7.520-win64-signed"
test ! -e "$CAP_DIR/modules/TRASH_2/dep/new.hor"
test ! -e "$CAP_DIR/modules/TRASH_2/temp"
appimage=$(find "$CAP_DIR" -type f -name '*.AppImage' -print -quit)
if [ -n "$appimage" ]; then
    echo "AppImage artifact must not be present in the OCI image: $appimage" >&2
    exit 1
fi

python3 - <<'PY'
import sklearn
import xgboost
assert sklearn.__version__ == "1.7.2"
assert xgboost.__version__ == "3.0.5"
PY
Rscript -e 'library(Biostrings); library(GenomicRanges); library(msa)'

rm -rf "$RESULTS_DIR" "$WORK_DIR"
/usr/local/bin/cap \
    --assembly "$CAP_DIR/test/ath_Chr1_extraction_trc.fasta" \
    --cores 1 \
    --outdir "$RESULTS_DIR" \
    -work-dir "$WORK_DIR"

file_count=$(find "$RESULTS_DIR" -maxdepth 1 -type f | wc -l)
test "$file_count" -eq 16 || {
    echo "Expected 16 outputs, found $file_count." >&2
    exit 1
}

# Rasterized plots can vary with the base OS font/rendering stack. Require the
# PNG to exist and be nonempty, but keep exact regression hashes for every
# deterministic data/model output.
find "$RESULTS_DIR" -maxdepth 1 -type f -name '*.png' -size +0c \
    -print -quit | grep -q .

grep -v '[.]png$' "$EXPECTED" > /tmp/expected-non-png.sha256
(
    cd "$RESULTS_DIR"
    sha256sum --check /tmp/expected-non-png.sha256
)
rm /tmp/expected-non-png.sha256

echo "OCI image test passed: 16 outputs and all 15 non-PNG hashes match."
