#!/usr/bin/env bash
#SBATCH --job-name=CAP
#SBATCH --cpus-per-task=1
#SBATCH --mem=36G
#SBATCH --output=CAP-%j.out
#SBATCH --error=CAP-%j.err

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  sbatch scripts/run-cap-slurm.sh --assembly FILE --outdir DIR [options]

Required:
  --assembly FILE       Genome assembly in FASTA format
  --outdir DIR          Results directory on persistent/shared storage

Optional:
  --trash2 DIR          Valid precomputed TRASH_2 output directory
  --te-gff FILE         EDTA transposable-element annotation
  --gene-gff FILE       Helixer gene annotation
  --metadata FILE       Existing chromosome metadata CSV
  --templates FILE      TRASH_2 repeat templates
  --max-rep-size INT    Maximum tandem-repeat size
  --min-rep-size INT    Minimum tandem-repeat size
  --env-name NAME       Conda environment (default: cap-pipeline)
  --conda EXE           Conda executable (default: conda from PATH)
  --help                Show this help

Run outside SLURM only for explicit local validation:
  CAP_SLURM_LOCAL_TEST=1 bash scripts/run-cap-slurm.sh ...
EOF
}

fail() {
    echo "Error: $*" >&2
    exit 1
}

require_value() {
    [ "$#" -ge 2 ] && [ -n "$2" ] || fail "$1 requires a value."
}

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
ENV_NAME="${CAP_ENV_NAME:-cap-pipeline}"
CONDA_EXE="${CAP_CONDA_EXE:-conda}"
ASSEMBLY=""
OUTDIR=""
TRASH2=""
TE_GFF=""
GENE_GFF=""
METADATA=""
TEMPLATES=""
MAX_REP_SIZE=""
MIN_REP_SIZE=""

while [ "$#" -gt 0 ]; do
    case "$1" in
        --assembly)
            require_value "$@"; ASSEMBLY=$2; shift 2 ;;
        --outdir)
            require_value "$@"; OUTDIR=$2; shift 2 ;;
        --trash2)
            require_value "$@"; TRASH2=$2; shift 2 ;;
        --te-gff)
            require_value "$@"; TE_GFF=$2; shift 2 ;;
        --gene-gff)
            require_value "$@"; GENE_GFF=$2; shift 2 ;;
        --metadata)
            require_value "$@"; METADATA=$2; shift 2 ;;
        --templates)
            require_value "$@"; TEMPLATES=$2; shift 2 ;;
        --max-rep-size)
            require_value "$@"; MAX_REP_SIZE=$2; shift 2 ;;
        --min-rep-size)
            require_value "$@"; MIN_REP_SIZE=$2; shift 2 ;;
        --env-name)
            require_value "$@"; ENV_NAME=$2; shift 2 ;;
        --conda)
            require_value "$@"; CONDA_EXE=$2; shift 2 ;;
        --help|-h)
            usage; exit 0 ;;
        *)
            fail "unknown argument: $1" ;;
    esac
done

[ -n "$ASSEMBLY" ] || fail "--assembly is required."
[ -n "$OUTDIR" ] || fail "--outdir is required."
[ -n "${SLURM_JOB_ID:-}" ] || [ "${CAP_SLURM_LOCAL_TEST:-0}" = 1 ] || \
    fail "submit this script with sbatch, or set CAP_SLURM_LOCAL_TEST=1 for local validation."

[ -f "$ASSEMBLY" ] || fail "assembly FASTA not found: $ASSEMBLY"
[ -z "$TE_GFF" ] || [ -f "$TE_GFF" ] || fail \
    "te annotation not found: $TE_GFF"
[ -z "$GENE_GFF" ] || [ -f "$GENE_GFF" ] || fail \
    "gene annotation not found: $GENE_GFF"
[ -z "$METADATA" ] || [ -f "$METADATA" ] || fail \
    "metadata CSV not found: $METADATA"
[ -z "$TEMPLATES" ] || [ -f "$TEMPLATES" ] || fail \
    "template file not found: $TEMPLATES"

ASSEMBLY="$(cd -- "$(dirname -- "$ASSEMBLY")" && pwd)/$(basename -- "$ASSEMBLY")"
OUTDIR_PARENT="$(dirname -- "$OUTDIR")"
mkdir -p -- "$OUTDIR_PARENT"
OUTDIR="$(cd -- "$OUTDIR_PARENT" && pwd)/$(basename -- "$OUTDIR")"
[ -z "$TE_GFF" ] || TE_GFF="$(cd -- "$(dirname -- "$TE_GFF")" && pwd)/$(basename -- "$TE_GFF")"
[ -z "$GENE_GFF" ] || GENE_GFF="$(cd -- "$(dirname -- "$GENE_GFF")" && pwd)/$(basename -- "$GENE_GFF")"
[ -z "$METADATA" ] || METADATA="$(cd -- "$(dirname -- "$METADATA")" && pwd)/$(basename -- "$METADATA")"
[ -z "$TEMPLATES" ] || TEMPLATES="$(cd -- "$(dirname -- "$TEMPLATES")" && pwd)/$(basename -- "$TEMPLATES")"

for item in "max-rep-size:$MAX_REP_SIZE" "min-rep-size:$MIN_REP_SIZE"; do
    label=${item%%:*}
    value=${item#*:}
    case "$value" in
        "") ;;
        *[!0-9]*|0) fail "--$label must be a positive integer." ;;
    esac
done
if [ -n "$MIN_REP_SIZE" ] && [ -n "$MAX_REP_SIZE" ] && \
    [ "$MIN_REP_SIZE" -gt "$MAX_REP_SIZE" ]; then
    fail "--min-rep-size must not exceed --max-rep-size."
fi

if [ -n "$TRASH2" ]; then
    [ -d "$TRASH2" ] || fail "TRASH_2 directory not found: $TRASH2"
    TRASH2="$(cd -- "$TRASH2" && pwd)"
    assembly_name=$(basename -- "$ASSEMBLY")
    repeats="$TRASH2/${assembly_name}_repeats_with_seq.csv"
    arrays="$TRASH2/${assembly_name}_arrays.csv"
    [ -f "$repeats" ] && [ -f "$arrays" ] || fail \
        "invalid --trash2 directory; expected $repeats and $arrays"
fi

if [[ "$CONDA_EXE" == */* ]]; then
    [ -x "$CONDA_EXE" ] || fail "Conda executable is not executable: $CONDA_EXE"
else
    command -v "$CONDA_EXE" >/dev/null 2>&1 || fail \
        "Conda executable not found in PATH: $CONDA_EXE"
fi

CORES="${SLURM_CPUS_PER_TASK:-1}"
case "$CORES" in
    *[!0-9]*|0) fail "SLURM_CPUS_PER_TASK must be a positive integer." ;;
esac

JOB_ID="${SLURM_JOB_ID:-local-$$}"
TMP_BASE="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
[ -d "$TMP_BASE" ] || fail "temporary directory base not found: $TMP_BASE"
WORK_DIR=$(mktemp -d "$TMP_BASE/cap-nextflow-$JOB_ID-XXXXXXXX")

mkdir -p -- "$OUTDIR"

nextflow_args=(
    nextflow run "$PROJECT_DIR"
    --assembly "$ASSEMBLY"
    --cores "$CORES"
    --outdir "$OUTDIR"
    -work-dir "$WORK_DIR"
)
[ -z "$TRASH2" ] || nextflow_args+=(--trash2 "$TRASH2")
[ -z "$TE_GFF" ] || nextflow_args+=(--te_gff "$TE_GFF")
[ -z "$GENE_GFF" ] || nextflow_args+=(--gene_gff "$GENE_GFF")
[ -z "$METADATA" ] || nextflow_args+=(--metadata "$METADATA")
[ -z "$TEMPLATES" ] || nextflow_args+=(--templates "$TEMPLATES")
[ -z "$MAX_REP_SIZE" ] || nextflow_args+=(--max_rep_size "$MAX_REP_SIZE")
[ -z "$MIN_REP_SIZE" ] || nextflow_args+=(--min_rep_size "$MIN_REP_SIZE")

command=("$CONDA_EXE" run -n "$ENV_NAME" "${nextflow_args[@]}")

printf 'CAP job ID: %s\n' "$JOB_ID"
printf 'CPUs: %s\n' "$CORES"
printf 'Work directory: %s\n' "$WORK_DIR"
printf 'Command:'
printf ' %q' "${command[@]}"
printf '\n'

if [ "${CAP_SLURM_DRY_RUN:-0}" = 1 ]; then
    status=0
else
    status=0
    "${command[@]}" || status=$?
fi

if [ "$status" -eq 0 ]; then
    printf 'CAP completed successfully.\n'
else
    printf 'CAP failed; retain the work directory for diagnosis.\n' >&2
fi
printf 'Nextflow work directory retained at: %s\n' "$WORK_DIR" >&2
printf 'Remove it when no longer needed with:\n  rm -rf -- %q\n' "$WORK_DIR" >&2
exit "$status"
