#!/usr/bin/env bash
#SBATCH --job-name=CAP-AppImage
#SBATCH --cpus-per-task=1
#SBATCH --mem=36G
#SBATCH --output=CAP-AppImage-%j.out
#SBATCH --error=CAP-AppImage-%j.err

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  sbatch scripts/run-cap-appimage-slurm.sh --appimage FILE --assembly FILE --outdir DIR [options]

Required:
  --appimage FILE       Executable CAP AppImage or extracted AppRun
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
  --help                Show this help

Run outside SLURM only for explicit local validation:
  CAP_SLURM_LOCAL_TEST=1 bash scripts/run-cap-appimage-slurm.sh ...
EOF
}

fail() {
    echo "Error: $*" >&2
    exit 1
}

require_value() {
    [ "$#" -ge 2 ] && [ -n "$2" ] || fail "$1 requires a value."
}

canonical_file() {
    local value=$1
    local label=$2
    [ -f "$value" ] || fail "$label not found: $value"
    printf '%s/%s\n' "$(cd -- "$(dirname -- "$value")" && pwd)" \
        "$(basename -- "$value")"
}

APPIMAGE="${CAP_APPIMAGE:-}"
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
        --appimage)
            require_value "$@"; APPIMAGE=$2; shift 2 ;;
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
        --help|-h)
            usage; exit 0 ;;
        *)
            fail "unknown argument: $1" ;;
    esac
done

[ -n "$APPIMAGE" ] || fail "--appimage is required."
[ -n "$ASSEMBLY" ] || fail "--assembly is required."
[ -n "$OUTDIR" ] || fail "--outdir is required."
[ -n "${SLURM_JOB_ID:-}" ] || [ "${CAP_SLURM_LOCAL_TEST:-0}" = 1 ] || \
    fail "submit this script with sbatch, or set CAP_SLURM_LOCAL_TEST=1 for local validation."

APPIMAGE=$(canonical_file "$APPIMAGE" "CAP AppImage or AppRun")
[ -x "$APPIMAGE" ] || fail "CAP AppImage or AppRun is not executable: $APPIMAGE"
ASSEMBLY=$(canonical_file "$ASSEMBLY" "assembly FASTA")
[ -z "$TE_GFF" ] || TE_GFF=$(canonical_file "$TE_GFF" "te annotation")
[ -z "$GENE_GFF" ] || GENE_GFF=$(canonical_file "$GENE_GFF" "gene annotation")
[ -z "$METADATA" ] || METADATA=$(canonical_file "$METADATA" "metadata CSV")
[ -z "$TEMPLATES" ] || TEMPLATES=$(canonical_file "$TEMPLATES" "template file")

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
    TRASH2=$(cd -- "$TRASH2" && pwd)
    assembly_name=$(basename -- "$ASSEMBLY")
    repeats="$TRASH2/${assembly_name}_repeats_with_seq.csv"
    arrays="$TRASH2/${assembly_name}_arrays.csv"
    [ -f "$repeats" ] && [ -f "$arrays" ] || fail \
        "invalid --trash2 directory; expected $repeats and $arrays"
fi

CORES=${SLURM_CPUS_PER_TASK:-1}
case "$CORES" in
    *[!0-9]*|0) fail "SLURM_CPUS_PER_TASK must be a positive integer." ;;
esac

OUTDIR_PARENT=$(dirname -- "$OUTDIR")
mkdir -p -- "$OUTDIR_PARENT"
OUTDIR="$(cd -- "$OUTDIR_PARENT" && pwd)/$(basename -- "$OUTDIR")"

JOB_ID=${SLURM_JOB_ID:-local-$$}
TMP_BASE=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
[ -d "$TMP_BASE" ] || fail "temporary directory base not found: $TMP_BASE"
WORK_DIR=$(mktemp -d "$TMP_BASE/cap-appimage-nextflow-$JOB_ID-XXXXXXXX")

mkdir -p -- "$OUTDIR"
command=(
    "$APPIMAGE"
    --assembly "$ASSEMBLY"
    --cores "$CORES"
    --outdir "$OUTDIR"
    -work-dir "$WORK_DIR"
)
[ -z "$TRASH2" ] || command+=(--trash2 "$TRASH2")
[ -z "$TE_GFF" ] || command+=(--te_gff "$TE_GFF")
[ -z "$GENE_GFF" ] || command+=(--gene_gff "$GENE_GFF")
[ -z "$METADATA" ] || command+=(--metadata "$METADATA")
[ -z "$TEMPLATES" ] || command+=(--templates "$TEMPLATES")
[ -z "$MAX_REP_SIZE" ] || command+=(--max_rep_size "$MAX_REP_SIZE")
[ -z "$MIN_REP_SIZE" ] || command+=(--min_rep_size "$MIN_REP_SIZE")

printf 'CAP AppImage job ID: %s\n' "$JOB_ID"
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
    printf 'CAP AppImage completed successfully.\n'
else
    printf 'CAP AppImage failed; retain the work directory for diagnosis.\n' >&2
fi
printf 'Nextflow work directory retained at: %s\n' "$WORK_DIR" >&2
printf 'Remove it when no longer needed with:\n  rm -rf -- %q\n' "$WORK_DIR" >&2
exit "$status"
