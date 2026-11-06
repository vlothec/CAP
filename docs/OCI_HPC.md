# Running CAP containers on HPC systems

CAP's supported scheduler architecture uses one SLURM allocation for the whole
workflow. Nextflow keeps its local executor inside that allocation and must not
submit individual pipeline processes to SLURM.

```text
sbatch
└── one allocation on one node
    └── Docker, Podman, or Apptainer
        └── CAP image
            └── Nextflow local executor
```

Partition names, accounts, time limits, and filesystem paths are site-specific.
Supply them as `sbatch` options rather than hard-coding them in CAP.

## Before choosing a runtime

Ask the cluster administrators which container runtime is supported on compute
nodes. Docker is commonly prohibited because its daemon is privileged. Rootless
Podman may be available, while Apptainer/Singularity is the most common HPC
choice. Also determine whether compute nodes can reach GHCR; if not, pull or
convert the image in advance on a connected system.

The convenient release reference is
`ghcr.io/mellifluousjdg/cap:1.0.3`. For an exact immutable image, use
`ghcr.io/mellifluousjdg/cap@sha256:2563036a264a3d7b72bba57215229df599af15a36ae0c4d7155b70a3127b6f95`.

## Rootless Podman

Create a final-results directory on durable shared storage, and use node-local
storage for temporary Nextflow work:

```bash
#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=36G

set -euo pipefail

image=ghcr.io/mellifluousjdg/cap:1.0.3
assembly=/shared/data/genome.fasta
outdir=/shared/results/genome-cap
tmp_base="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
job_id=${SLURM_JOB_ID:?submit this script with sbatch}
work_root=$(mktemp -d "$tmp_base/cap-$job_id-XXXXXXXX")

mkdir -p "$outdir"

podman run --rm --userns=keep-id \
  -v "$(dirname "$assembly"):/data:ro" \
  -v "$outdir:/results" \
  -v "$work_root:/work" \
  "$image" \
  --assembly "/data/$(basename "$assembly")" \
  --outdir /results \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  -work-dir /work/nextflow

printf 'Work directory retained at: %s\n' "$work_root"
printf 'Remove it when no longer needed with:\n  rm -rf -- %q\n' "$work_root"
```

`--userns=keep-id` is important for rootless Podman: it makes bind-mounted
outputs belong to the submitting host user instead of a subordinate UID.
SELinux-enabled hosts may require `:Z` on private bind mounts; consult the
cluster's container-runtime policy before relabeling shared paths.

## Docker

Use Docker only where administrators explicitly support it on compute nodes.
The command shape is similar, but map the runtime process to the submitting
user so bind-mounted outputs remain host-owned:

```bash
work_root="${SLURM_TMPDIR}/cap-work"
mkdir -p "$work_root/home" /shared/results/genome-cap

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -e HOME=/work/home \
  -e NXF_HOME=/work/home/.nextflow \
  -v /shared/data:/data:ro \
  -v /shared/results/genome-cap:/results \
  -v "$work_root:/work" \
  ghcr.io/mellifluousjdg/cap:1.0.3 \
  --assembly /data/genome.fasta \
  --outdir /results \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  -work-dir /work/nextflow
```

The example creates the writable home and result directories before launching.
Its work directory is retained; keep it for diagnosis or resume, and remove it
manually when it is no longer needed.

## Apptainer or Singularity

Because the image is public, a connected login or transfer node can normally
convert it to a SIF file:

```bash
apptainer pull cap-1.0.3.sif \
  docker://ghcr.io/mellifluousjdg/cap:1.0.3
```

Then copy the immutable SIF to shared storage and execute it in one allocation:

```bash
outdir=/shared/results/genome-cap
tmp_base="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
job_id=${SLURM_JOB_ID:?run this command inside a SLURM allocation}
work_root=$(mktemp -d "$tmp_base/cap-$job_id-XXXXXXXX")

mkdir -p "$outdir" "$work_root/home"

apptainer exec \
  --env HOME=/work/home \
  --env NXF_HOME=/work/home/.nextflow \
  --bind /shared/data:/data:ro \
  --bind "$outdir:/results" \
  --bind "$work_root:/work" \
  cap-1.0.3.sif \
  /usr/local/bin/cap \
  --assembly /data/genome.fasta \
  --outdir /results \
  --cores "${SLURM_CPUS_PER_TASK:-1}" \
  -work-dir /work/nextflow

printf 'Work directory retained at: %s\n' "$work_root"
printf 'Remove it when no longer needed with:\n  rm -rf -- %q\n' "$work_root"
```

Apptainer does not necessarily honor an OCI image's entrypoint, so invoke
`/usr/local/bin/cap` explicitly. Final results remain on shared storage. The
node-local work directory is retained for diagnosis or Nextflow resume; remove
it manually when it is no longer needed.

## Offline transfer

For clusters without registry access, export an OCI archive on a connected
host:

```bash
podman save --format oci-archive \
  -o cap-1.0.3.oci.tar \
  ghcr.io/mellifluousjdg/cap:1.0.3
sha256sum cap-1.0.3.oci.tar > cap-1.0.3.oci.tar.sha256
```

Transfer both files, verify the checksum, and import with the cluster's
supported runtime. An OCI archive is not a substitute for an Apptainer SIF when
the cluster only supports Apptainer.
