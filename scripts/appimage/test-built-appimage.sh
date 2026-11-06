#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE=${1:-}
[ -n "$IMAGE" ] || {
    echo "Usage: $0 path/to/CAP-x86_64.AppImage" >&2
    exit 1
}
IMAGE="$(cd -- "$(dirname -- "$IMAGE")" && pwd)/$(basename -- "$IMAGE")"
[ -x "$IMAGE" ] || {
    echo "AppImage is not executable: $IMAGE" >&2
    exit 1
}

TEST_DIR="${APPIMAGE_TEST_DIR:-$PROJECT_DIR/build/appimage/artifact-test}"
EXPECTED="$TEST_DIR/expected-non-png.sha256"
ASSEMBLY="$PROJECT_DIR/test/ath_Chr1_extraction_trc.fasta"
rm -rf "$TEST_DIR"
mkdir -p "$TEST_DIR/direct/cwd" "$TEST_DIR/direct/home" \
    "$TEST_DIR/extracted/cwd" "$TEST_DIR/extracted/home" \
    "$TEST_DIR/unpack"
grep -v '[.]png$' "$PROJECT_DIR/test/expected-results.sha256" > "$EXPECTED"

verify_results() {
    local mode=$1
    local results=$2
    local file_count

    file_count=$(find "$results" -maxdepth 1 -type f | wc -l)
    test "$file_count" -eq 16 || {
        echo "$mode: expected 16 outputs, found $file_count." >&2
        exit 1
    }
    find "$results" -maxdepth 1 -type f -name '*.png' -size +0c \
        -print -quit | grep -q .
    (cd "$results" && sha256sum --check "$EXPECTED")
    echo "$mode passed: 16 outputs and all 15 non-PNG hashes match."
}

# Exercise the AppImage's normal FUSE-mounted execution path. Do not permit the
# AppImage runtime's extraction fallback in this part of the test.
(
    unset APPIMAGE_EXTRACT_AND_RUN APPIMAGE_EXTRACT_AND_RUN_ONCE APPDIR
    export HOME="$TEST_DIR/direct/home"
    cd "$TEST_DIR/direct/cwd"
    if ! "$IMAGE" \
        --assembly "$ASSEMBLY" \
        --cores 1 \
        --outdir "$TEST_DIR/direct/results" \
        -work-dir "$TEST_DIR/direct/work"; then
        echo "Normal AppImage execution failed. This test requires working FUSE 2 support (including libfuse.so.2) on the host." >&2
        exit 1
    fi
)
verify_results "Normal AppImage execution" "$TEST_DIR/direct/results"

RUNTIME="$TEST_DIR/extracted/runtime"
(
    cd "$TEST_DIR/unpack"
    "$IMAGE" --appimage-extract >/dev/null
    mv squashfs-root "$RUNTIME"
)

export APPDIR="$RUNTIME"
export HOME="$TEST_DIR/extracted/home"
"$RUNTIME/AppRun" --shell -c '
    set -e
    for command in nextflow java Rscript python3 mafft nhmmer; do
        command -v "$command" >/dev/null
    done
    test "$(Rscript -e "cat(R.home())")" = "$APPDIR/usr/lib/R"
    Rscript -e "library(Biostrings); library(GenomicRanges); library(msa)"
    python3 -c "import sklearn, xgboost; assert sklearn.__version__ == '\''1.7.2'\''; assert xgboost.__version__ == '\''3.0.5'\''"
    test -x "$APPDIR/usr/lib/cap/bin/ctw-calc"
    ldd "$APPDIR/usr/lib/cap/bin/ctw-calc" | grep -F "$APPDIR/usr/lib/libstdc++.so.6"
    ldd "$APPDIR/usr/lib/cap/bin/ctw-calc" | grep -F "$APPDIR/usr/lib/libgcc_s.so.1"
    "$APPDIR/usr/lib/cap/bin/ctw-calc" >/dev/null 2>&1 || test "$?" -eq 1
    test ! -e "$APPDIR/usr/lib/cap/modules/TRASH_2/dep/hmmer"
    test ! -e "$APPDIR/usr/lib/cap/modules/TRASH_2/dep/mafft-7.520-win64-signed"
'

cd "$TEST_DIR/extracted/cwd"
"$RUNTIME/AppRun" \
    --assembly "$ASSEMBLY" \
    --cores 1 \
    --outdir "$TEST_DIR/extracted/results" \
    -work-dir "$TEST_DIR/extracted/work"
verify_results "Extracted AppImage execution" "$TEST_DIR/extracted/results"

echo "AppImage artifact test passed in normal and extracted execution modes."
