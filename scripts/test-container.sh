#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_CACHE=$(mktemp -d)
trap 'rm -rf "$PYTHON_CACHE"' EXIT
export PYTHONPYCACHEPREFIX="$PYTHON_CACHE"
bash "$PROJECT_DIR/scripts/compliance/test-source-obligations.sh"
python3 "$PROJECT_DIR/scripts/compliance/test-source-review-policy.py"
bash "$PROJECT_DIR/scripts/compliance/test-collect-conda-sources.sh"
bash "$PROJECT_DIR/scripts/compliance/test-collect-conda-notices.sh"
python3 "$PROJECT_DIR/scripts/compliance/test-merge-conda-source-corpora.py"
python3 "$PROJECT_DIR/scripts/compliance/test-collect-debian-sources.py"
bash "$PROJECT_DIR/scripts/compliance/test-conda-source-notices.sh"
bash "$PROJECT_DIR/scripts/compliance/test-finalize-compliance-checksums.sh"
bash "$PROJECT_DIR/scripts/compliance/test-resolve-release-revisions.sh"
bash "$PROJECT_DIR/scripts/compliance/test-sbom.sh"

runtime_lock="$PROJECT_DIR/conda-runtime-linux-64.lock"
full_lock="$PROJECT_DIR/conda-linux-64.lock"
test -s "$runtime_lock"
while IFS= read -r package_url; do
    grep -Fqx "$package_url" "$full_lock" || {
        echo "Runtime lock URL is absent from the validated full lock: $package_url" >&2
        exit 1
    }
done < <(grep '^https://' "$runtime_lock")
for build_package in gcc_impl_linux-64 gxx_impl_linux-64 \
    gfortran_impl_linux-64 binutils_impl_linux-64 sysroot_linux-64 \
    kernel-headers_linux-64 cmake make git; do
    if grep -Eq "/${build_package}-[0-9]" "$runtime_lock"; then
        echo "Build package present in runtime lock: $build_package" >&2
        exit 1
    fi
done

test -f "$PROJECT_DIR/modules/TRASH_2/src/TRASH.R" || {
    echo "TRASH_2 submodule is not initialized." >&2
    exit 1
}
for license_file in \
    "$PROJECT_DIR/THIRD_PARTY_NOTICES.md" \
    "$PROJECT_DIR/docs/DEPENDENCY_SOURCE_REVIEW.md" \
    "$PROJECT_DIR/bin/src/BCT/README.md" \
    "$PROJECT_DIR/licenses/GPL-2.0.txt" \
    "$PROJECT_DIR/licenses/GPL-3.0.txt"; do
    test -s "$license_file" || {
        echo "Required licensing file missing or empty: $license_file" >&2
        exit 1
    }
done

echo "OCI container static checks passed."
