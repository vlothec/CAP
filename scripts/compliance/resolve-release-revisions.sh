#!/usr/bin/env bash

set -euo pipefail

image_cap_revision=${1-}
image_trash_revision=${2-}
checkout_cap_revision=${3-}
checkout_trash_revision=${4-}
expected_cap_revision=${CAP_EXPECTED_REVISION-}
expected_trash_revision=${CAP_EXPECTED_TRASH2_REVISION-}

if [[ -n "$expected_cap_revision" || -n "$expected_trash_revision" ]]; then
    [[ -n "$expected_cap_revision" && -n "$expected_trash_revision" ]] || {
        echo "Both explicit expected revisions must be provided together." >&2
        exit 1
    }
    [[ "$checkout_cap_revision" == "$expected_cap_revision" ]] || {
        echo "Checkout revision does not match explicit CAP revision." >&2
        exit 1
    }
    [[ "$checkout_trash_revision" == "$expected_trash_revision" ]] || {
        echo "Checkout TRASH_2 revision does not match explicit revision." >&2
        exit 1
    }
    printf '%s\n%s\n' "$expected_cap_revision" "$expected_trash_revision"
    exit 0
fi

if [[ -z "$image_cap_revision" || "$image_cap_revision" == "UNSPECIFIED" ]]; then
    echo "Image has no specific CAP revision label." >&2
    exit 1
fi
if [[ -z "$image_trash_revision" || \
      "$image_trash_revision" == "UNSPECIFIED" ]]; then
    echo "Image has no specific TRASH_2 revision label." >&2
    exit 1
fi
[[ "$checkout_cap_revision" == "$image_cap_revision" ]] || {
    echo "Checkout revision does not match image CAP revision." >&2
    exit 1
}
[[ "$checkout_trash_revision" == "$image_trash_revision" ]] || {
    echo "Checkout TRASH_2 revision does not match image revision." >&2
    exit 1
}
printf '%s\n%s\n' "$image_cap_revision" "$image_trash_revision"
