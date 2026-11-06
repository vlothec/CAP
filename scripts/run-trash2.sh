#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/run-trash2.sh --assembly FILE --outdir DIR [options]

Run CAP's pinned TRASH_2 revision in an existing CAP Conda environment.

Required:
  --assembly FILE       Input FASTA file
  --outdir DIR          Output directory; created if its parent exists

Options:
  --env-name NAME       Conda environment name (default: cap-pipeline)
  --conda EXE           Conda executable (default: conda from PATH)
  --cores N             Worker count (default: 1)
  --max-rep-size N      Maximum repeat size (default: 1000)
  --min-rep-size N      Minimum repeat size (default: 7)
  --templates FILE      Optional repeat-template FASTA
  -h, --help            Show this help
EOF
}

fail() {
    echo "run-trash2: $*" >&2
    exit 2
}

require_value() {
    [ "$#" -ge 2 ] && [ -n "$2" ] || fail "option $1 requires a value"
}

require_positive_integer() {
    case "$2" in
        ''|*[!0-9]*|0) fail "$1 must be a positive integer: $2" ;;
    esac
}

canonical_existing_file() {
    local value=$1
    local label=$2
    [ -f "$value" ] || fail "$label file not found: $value"
    readlink -f -- "$value"
}

canonical_output_dir() {
    local value=$1
    local parent
    local name

    [ -n "$value" ] || fail "output directory must not be empty"
    if [ -e "$value" ] && [ ! -d "$value" ]; then
        fail "output path exists but is not a directory: $value"
    fi
    if [ -d "$value" ]; then
        readlink -f -- "$value"
        return
    fi

    parent=$(dirname -- "$value")
    name=$(basename -- "$value")
    [ "$name" != "." ] && [ "$name" != ".." ] ||
        fail "invalid output directory: $value"
    [ -d "$parent" ] || fail "output parent directory not found: $parent"
    printf '%s/%s\n' "$(readlink -f -- "$parent")" "$name"
}

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)
TRASH2_SCRIPT="$PROJECT_DIR/modules/TRASH_2/src/TRASH.R"
ENV_NAME=${CAP_ENV_NAME:-cap-pipeline}
CONDA_EXE=${CAP_CONDA_EXE:-conda}
ASSEMBLY=
OUTDIR=
CORES=1
MAX_REP_SIZE=1000
MIN_REP_SIZE=7
TEMPLATES=

while [ "$#" -gt 0 ]; do
    case "$1" in
        --assembly)
            require_value "$@"; ASSEMBLY=$2; shift 2 ;;
        --outdir)
            require_value "$@"; OUTDIR=$2; shift 2 ;;
        --env-name)
            require_value "$@"; ENV_NAME=$2; shift 2 ;;
        --conda)
            require_value "$@"; CONDA_EXE=$2; shift 2 ;;
        --cores)
            require_value "$@"; CORES=$2; shift 2 ;;
        --max-rep-size)
            require_value "$@"; MAX_REP_SIZE=$2; shift 2 ;;
        --min-rep-size)
            require_value "$@"; MIN_REP_SIZE=$2; shift 2 ;;
        --templates)
            require_value "$@"; TEMPLATES=$2; shift 2 ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift
            [ "$#" -eq 0 ] || fail "unexpected positional argument: $1"
            ;;
        -*) fail "unknown option: $1" ;;
        *) fail "unexpected positional argument: $1" ;;
    esac
done

[ -n "$ASSEMBLY" ] || fail "--assembly is required"
[ -n "$OUTDIR" ] || fail "--outdir is required"
[ -n "$ENV_NAME" ] || fail "environment name must not be empty"
require_positive_integer --cores "$CORES"
require_positive_integer --max-rep-size "$MAX_REP_SIZE"
require_positive_integer --min-rep-size "$MIN_REP_SIZE"
[ "$MIN_REP_SIZE" -le "$MAX_REP_SIZE" ] ||
    fail "--min-rep-size must not exceed --max-rep-size"
[ -f "$TRASH2_SCRIPT" ] ||
    fail "TRASH_2 submodule is not initialized; run git submodule update --init --recursive"
if [[ "$CONDA_EXE" == */* ]]; then
    [ -x "$CONDA_EXE" ] || fail "Conda executable is not executable: $CONDA_EXE"
else
    command -v "$CONDA_EXE" >/dev/null 2>&1 ||
        fail "Conda executable not found in PATH: $CONDA_EXE"
fi

ASSEMBLY=$(canonical_existing_file "$ASSEMBLY" assembly)
OUTDIR=$(canonical_output_dir "$OUTDIR")
if [ -n "$TEMPLATES" ]; then
    TEMPLATES=$(canonical_existing_file "$TEMPLATES" templates)
fi

conda_environment=$(
    "$CONDA_EXE" env list | awk -v requested="$ENV_NAME" '
        NF && $1 !~ /^#/ {
            if ($1 == requested) {
                print requested
                exit
            }
            path = $NF
            count = split(path, parts, "/")
            if (parts[count] == requested) {
                print path
                exit
            }
        }
    '
)
[ -n "$conda_environment" ] ||
    fail "Conda environment '$ENV_NAME' does not exist; run make install ENV_NAME=$ENV_NAME"
if [[ "$conda_environment" == */* ]]; then
    conda_environment_args=(-p "$conda_environment")
else
    conda_environment_args=(-n "$conda_environment")
fi

required_r_packages='stringr,stringdist,seqinr,doParallel,getopt,ape,Biostrings,msa'
startup_file=$(mktemp)
trap 'rm -f "$startup_file"' EXIT
cat > "$startup_file" <<'RSCRIPT'
install.packages <- function(...) {
  stop("runtime package installation is disabled; recreate CAP with make install",
       call. = FALSE)
}
RSCRIPT

preflight=$(cat <<'RSCRIPT'
required <- strsplit(Sys.getenv("CAP_TRASH2_REQUIRED_R_PACKAGES"), ",", fixed = TRUE)[[1]]
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("missing locked R packages: ", paste(missing, collapse = ", "), call. = FALSE)
}
for (tool in c("mafft", "nhmmer")) {
  if (!nzchar(Sys.which(tool))) stop("missing locked executable: ", tool, call. = FALSE)
}
cat("Locked TRASH_2 dependencies are available.\n")
RSCRIPT
)

CAP_TRASH2_REQUIRED_R_PACKAGES="$required_r_packages" \
R_PROFILE_USER="$startup_file" \
http_proxy='http://127.0.0.1:9' \
https_proxy='http://127.0.0.1:9' \
HTTP_PROXY='http://127.0.0.1:9' \
HTTPS_PROXY='http://127.0.0.1:9' \
R_REPOSITORIES='@CRAN@=file:///nonexistent/cap-no-runtime-installs' \
    "$CONDA_EXE" run --no-capture-output "${conda_environment_args[@]}" \
    Rscript -e "$preflight"

mkdir -p -- "$OUTDIR"
arguments=(
    -f "$ASSEMBLY"
    -o "$OUTDIR"
    --cores_no "$CORES"
    --max_rep_size "$MAX_REP_SIZE"
    --min_rep_size "$MIN_REP_SIZE"
)
[ -z "$TEMPLATES" ] || arguments+=(--templates "$TEMPLATES")

R_PROFILE_USER="$startup_file" \
http_proxy='http://127.0.0.1:9' \
https_proxy='http://127.0.0.1:9' \
HTTP_PROXY='http://127.0.0.1:9' \
HTTPS_PROXY='http://127.0.0.1:9' \
R_REPOSITORIES='@CRAN@=file:///nonexistent/cap-no-runtime-installs' \
    "$CONDA_EXE" run --no-capture-output "${conda_environment_args[@]}" \
    Rscript "$TRASH2_SCRIPT" "${arguments[@]}"
