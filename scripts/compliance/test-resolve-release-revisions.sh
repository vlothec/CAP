#!/usr/bin/env bash

set -euo pipefail

PROJECT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
SCRIPT="$PROJECT_DIR/scripts/compliance/resolve-release-revisions.sh"
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

run_success() {
    local expected=$1
    shift
    "$@" > "$TMP_DIR/output"
    diff -u <(printf '%s' "$expected") "$TMP_DIR/output"
}

run_failure() {
    local message=$1
    shift
    set +e
    "$@" > "$TMP_DIR/output" 2> "$TMP_DIR/error"
    status=$?
    set -e
    [[ "$status" -ne 0 ]]
    grep -Fq "$message" "$TMP_DIR/error"
}

run_success $'cap-tag\ntrash-tag\n' \
    env CAP_EXPECTED_REVISION=cap-tag \
    CAP_EXPECTED_TRASH2_REVISION=trash-tag \
    "$SCRIPT" UNSPECIFIED UNSPECIFIED cap-tag trash-tag
run_failure 'explicit CAP revision' \
    env CAP_EXPECTED_REVISION=cap-tag \
    CAP_EXPECTED_TRASH2_REVISION=trash-tag \
    "$SCRIPT" UNSPECIFIED UNSPECIFIED wrong trash-tag
run_failure 'provided together' \
    env CAP_EXPECTED_REVISION=cap-tag \
    "$SCRIPT" cap-label trash-label cap-tag trash-tag
run_success $'cap-label\ntrash-label\n' \
    "$SCRIPT" cap-label trash-label cap-label trash-label
run_failure 'no specific TRASH_2 revision label' \
    "$SCRIPT" cap-label UNSPECIFIED cap-label trash-label
run_failure 'does not match image CAP revision' \
    "$SCRIPT" cap-label trash-label wrong trash-label

echo 'Compliance release revision resolution tests passed.'
