#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD_DIR="${APPIMAGE_BUILD_DIR:-$PROJECT_DIR/build/appimage}"
OUTPUT_DIR="${APPIMAGE_OUTPUT_DIR:-$PROJECT_DIR/dist}"
CACHE_DIR="${APPIMAGE_CACHE_DIR:-$BUILD_DIR/cache}"
APPDIR="$BUILD_DIR/CAP.AppDir"
CONDA_DIR="$BUILD_DIR/conda-bootstrap"
CTW_BUILD_DIR="$BUILD_DIR/ctw-build"

MINICONDA_VERSION=py313_25.7.0-2
MINICONDA_SHA256=dda3629462ba1cfa72eb74535214c2e315c77f1cfb0f02046537e99f1bf64abc
MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh"

# AppImageKit continuous asset retrieved 2026-10-04. The release URL is
# floating, but the accepted bytes are fixed by this required checksum.
APPIMAGETOOL_URL="${APPIMAGETOOL_URL:-https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage}"
APPIMAGETOOL_SHA256="${APPIMAGETOOL_SHA256:-b90f4a8b18967545fda78a445b27680a1642f1ef9488ced28b65398f2be7add2}"

fail() {
    echo "AppImage build error: $*" >&2
    exit 1
}

download_verified() {
    local url=$1
    local output=$2
    local checksum=$3

    if [ ! -f "$output" ]; then
        curl --fail --location --retry 3 "$url" --output "$output"
    fi
    echo "$checksum  $output" | sha256sum --check --strict
}

[ "$(uname -s)" = Linux ] || fail "only Linux builds are supported"
[ "$(uname -m)" = x86_64 ] || fail "only Linux x86-64 builds are supported"
for command in bash cmp curl git readelf sed sha256sum sort; do
    command -v "$command" >/dev/null || fail "required host command not found: $command"
done
[ -s "$PROJECT_DIR/conda-runtime-linux-64.lock" ] || fail "runtime lock is missing"
[ -s "$PROJECT_DIR/conda-linux-64.lock" ] || fail "development lock is missing"
[ -f "$PROJECT_DIR/modules/TRASH_2/src/TRASH.R" ] ||
    fail "TRASH_2 is not initialized; run git submodule update --init --recursive"
[ -n "$APPIMAGETOOL_SHA256" ] || fail \
    "set APPIMAGETOOL_SHA256 to the verified checksum for APPIMAGETOOL_URL"

mkdir -p "$CACHE_DIR" "$OUTPUT_DIR"
miniconda="$CACHE_DIR/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh"
appimagetool="$CACHE_DIR/appimagetool-x86_64.AppImage"
download_verified "$MINICONDA_URL" "$miniconda" "$MINICONDA_SHA256"
download_verified "$APPIMAGETOOL_URL" "$appimagetool" "$APPIMAGETOOL_SHA256"
chmod 0755 "$appimagetool"

rm -rf "$APPDIR" "$CONDA_DIR" "$CTW_BUILD_DIR"
bash "$miniconda" -b -p "$CONDA_DIR"
"$CONDA_DIR/bin/conda" create --yes --prefix "$APPDIR/usr" \
    --override-channels \
    --channel conda-forge \
    --channel bioconda \
    --file "$PROJECT_DIR/conda-runtime-linux-64.lock"

# Build CTW with the pinned development lock rather than the host compiler.
# Compile from a disposable source copy so the repository stays artifact-free.
"$CONDA_DIR/bin/conda" create --yes --prefix "$CTW_BUILD_DIR" \
    --override-channels \
    --channel conda-forge \
    --channel bioconda \
    --file "$PROJECT_DIR/conda-linux-64.lock"
mkdir -p "$CTW_BUILD_DIR/bin/src/BCT"
cp -a "$PROJECT_DIR/bin/src/BCT/." "$CTW_BUILD_DIR/bin/src/BCT/"
"$CONDA_DIR/bin/conda" run --prefix "$CTW_BUILD_DIR" \
    make -C "$CTW_BUILD_DIR/bin/src/BCT"
ctw_binary="$CTW_BUILD_DIR/bin/ctw-calc"
test -x "$ctw_binary" || fail "CTW build did not produce ctw-calc"

ctw_test_dir="$CTW_BUILD_DIR/smoke-test"
mkdir -p "$ctw_test_dir"
printf '>smoke\n' > "$ctw_test_dir/input.fa"
i=0
while [ "$i" -lt 30 ]; do
    printf 'ACGT' >> "$ctw_test_dir/input.fa"
    i=$((i + 1))
done
printf '\n' >> "$ctw_test_dir/input.fa"
"$ctw_binary" "$ctw_test_dir/input.fa" "$ctw_test_dir/output.csv" 100 \
    >/dev/null
printf 'chromosome,bin_mid,bin_value\nsmoke,61.00,20.57\n' > \
    "$ctw_test_dir/expected.csv"
cmp -s "$ctw_test_dir/expected.csv" "$ctw_test_dir/output.csv" ||
    fail "CTW execution smoke test failed"

highest_glibc=$(readelf --version-info "$ctw_binary" |
    grep -o 'GLIBC_[0-9.]*' | sort -V | tail -n 1)
[ "$highest_glibc" = GLIBC_2.14 ] ||
    fail "ctw-calc unexpectedly requires $highest_glibc"

mkdir -p "$APPDIR/usr/lib/cap" "$APPDIR/usr/share/applications" \
    "$APPDIR/usr/share/icons/hicolor/scalable/apps"
git -C "$PROJECT_DIR" archive --format=tar HEAD | \
    tar -xf - -C "$APPDIR/usr/lib/cap"
git -C "$PROJECT_DIR/modules/TRASH_2" archive --format=tar HEAD | \
    tar -xf - -C "$APPDIR/usr/lib/cap/modules/TRASH_2"
install -m 0755 "$ctw_binary" \
    "$APPDIR/usr/lib/cap/bin/ctw-calc"
rm -rf "$CTW_BUILD_DIR"
"$CONDA_DIR/bin/conda" clean --all --yes
rm -rf \
    "$APPDIR/usr/lib/cap/modules/TRASH_2/dep/hmmer" \
    "$APPDIR/usr/lib/cap/modules/TRASH_2/dep/mafft-7.520-win64-signed" \
    "$APPDIR/usr/lib/cap/modules/TRASH_2/dep/new.hor" \
    "$APPDIR/usr/lib/cap/modules/TRASH_2/temp"

install -m 0755 "$PROJECT_DIR/scripts/appimage/AppRun" "$APPDIR/AppRun"
install -m 0755 "$PROJECT_DIR/scripts/appimage/Rscript" "$APPDIR/usr/bin/Rscript"
install -m 0644 "$PROJECT_DIR/scripts/appimage/cap.desktop" \
    "$APPDIR/cap.desktop"
install -m 0644 "$PROJECT_DIR/scripts/appimage/cap.svg" "$APPDIR/cap.svg"
ln -s ../../../../../cap.svg \
    "$APPDIR/usr/share/icons/hicolor/scalable/apps/cap.svg"
install -m 0644 "$PROJECT_DIR/scripts/appimage/cap.desktop" \
    "$APPDIR/usr/share/applications/cap.desktop"

# The installed Rscript wrapper resolves R relative to itself. Conda's original
# Rscript embeds the build prefix and is deliberately replaced above.
test ! -e "$APPDIR/usr/bin/Rscript.real"
grep -F 'PREFIX="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"' \
    "$APPDIR/usr/bin/Rscript" >/dev/null

cap_revision="$(git -C "$PROJECT_DIR" rev-parse HEAD)"
trash_revision="$(git -C "$PROJECT_DIR/modules/TRASH_2" rev-parse HEAD)"
appimage_version="${APPIMAGE_VERSION:-${cap_revision:0:7}}"
if [[ ! "$appimage_version" =~ ^v?[0-9]+\.[0-9]+\.[0-9]+$ && \
      "$appimage_version" != "${cap_revision:0:7}" ]]; then
    fail "APPIMAGE_VERSION must be a semantic version such as v1.0.3"
fi
cat > "$APPDIR/usr/lib/cap/APPIMAGE-BUILD-METADATA.txt" <<EOF
CAP revision: $cap_revision
TRASH_2 revision: $trash_revision
Release version: $appimage_version
Runtime lock SHA-256: $(sha256sum "$PROJECT_DIR/conda-runtime-linux-64.lock" | cut -d' ' -f1)
Miniconda version: $MINICONDA_VERSION
AppImage architecture: x86_64
EOF

output="$OUTPUT_DIR/CAP-${appimage_version}-x86_64.AppImage"
rm -f "$output"
ARCH=x86_64 APPIMAGE_EXTRACT_AND_RUN=1 "$appimagetool" "$APPDIR" "$output"
chmod 0755 "$output"
sha256sum "$output" > "$output.sha256"
printf 'Built %s\n' "$output"
