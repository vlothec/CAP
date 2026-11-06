#!/usr/bin/env bash

set -euo pipefail

CAP_DIR=/opt/CAP

if [ "${1:-}" = "--shell" ]; then
    shift
    exec /bin/bash "$@"
fi

if [ "${1:-}" = "nextflow" ]; then
    exec "$@"
fi

exec nextflow run "$CAP_DIR" "$@"
