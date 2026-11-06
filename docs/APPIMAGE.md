# CAP AppImage

The CAP AppImage is a portable Linux x86-64 application containing CAP, its
pinned TRASH_2 module, Nextflow, Java, R/Bioconductor, Python, MAFFT, HMMER, the
prediction model, and CTW/BCT. It does not require Conda, Java, R, Python, or a
container runtime on the host.

## Download a release

Download these two assets from the CAP GitHub Release page:

```text
CAP-<version>-x86_64.AppImage
CAP-<version>-x86_64.AppImage.sha256
```

Keep both files in the same directory, verify the download, and make the
AppImage executable:

```bash
cd /path/to/downloads
sha256sum --check CAP-<version>-x86_64.AppImage.sha256
chmod +x CAP-<version>-x86_64.AppImage
```

Replace `<version>` with the version shown in the release. Use the checksum
published with that same release; do not reuse a checksum from another version.

The AppImage supports Linux x86-64 only. macOS, Windows, ARM64, and 32-bit Linux
are not supported by this artifact.

## Run CAP

Use absolute paths so input, output, and Nextflow work locations are
unambiguous:

```bash
appimage=/path/to/CAP-<version>-x86_64.AppImage
assembly=/path/to/genome.fasta
outdir=/path/to/cap-results
workdir=/path/to/cap-work

mkdir -p "$outdir" "$workdir"
"$appimage" \
  --assembly "$assembly" \
  --cores 1 \
  --outdir "$outdir" \
  -work-dir "$workdir"
```

Final CSV, RDS, model, and PNG files are written to `outdir`. Nextflow stores
intermediate files in `workdir`; remove that directory after a successful run
to reclaim space:

```bash
rm -rf /path/to/cap-work
```

Keep the work directory and `.nextflow.log` while diagnosing a failed run.

The number supplied to `--cores` should not exceed the CPUs available to the
process. CAP's default is one core. A larger assembly may also require more
memory and temporary storage.

### Optional inputs and parameters

The AppImage accepts the same workflow arguments as the Conda installation:

- `--te_gff`: EDTA transposable-element annotation in GFF3 format.
- `--gene_gff`: Helixer gene annotation in GFF format.
- `--metadata`: existing chromosome metadata CSV.
- `--templates`: repeat-template FASTA passed to TRASH_2.
- `--trash2`: existing TRASH_2 output directory.
- `--max_rep_size`: maximum tandem-repeat size; default `1000`.
- `--min_rep_size`: minimum tandem-repeat size; default `7`.

For example:

```bash
"$appimage" \
  --assembly /data/arabidopsis.fasta \
  --te_gff /data/EDTA.TEanno.split.gff3 \
  --gene_gff /data/helixer.gff \
  --cores 4 \
  --outdir /data/cap-results \
  -work-dir /data/cap-work
```

Run `"$appimage" --help` for Nextflow and pipeline help. See
`docs/PIPELINE_FLOW.md` for the complete output inventory and workflow details.

## Run on a system without FUSE

Normal execution may require FUSE support from the operating system. If direct
execution reports a FUSE, mount, or AppImage runtime error, extract the AppImage
and run the extracted application:

```bash
mkdir cap-appimage
cd cap-appimage
/path/to/CAP-<version>-x86_64.AppImage --appimage-extract

./squashfs-root/AppRun \
  --assembly /path/to/genome.fasta \
  --cores 1 \
  --outdir /path/to/cap-results \
  -work-dir /path/to/cap-work
```

Move the entire `squashfs-root` directory if desired; do not move individual
files out of it. CAP resolves its bundled runtime relative to that directory,
so extraction can be placed in any user-writable location. Remove the extracted
directory when it is no longer needed:

```bash
cd ..
rm -rf cap-appimage
```

## Run in one SLURM allocation

CAP includes a generic AppImage `sbatch` wrapper that starts the workflow once
inside one allocation. Nextflow uses its local executor and does not submit
separate jobs for pipeline processes.

Place the AppImage, assembly, and final output directory on storage visible to
the compute node, then submit:

```bash
sbatch scripts/run-cap-appimage-slurm.sh \
  --appimage /shared/apps/CAP-<version>-x86_64.AppImage \
  --assembly /shared/data/genome.fasta \
  --outdir /shared/results/genome-cap
```

The generic header requests one CPU and 36 GB of memory. Override resources with
site-appropriate `sbatch` options; the wrapper passes `SLURM_CPUS_PER_TASK` to
CAP:

```bash
sbatch --cpus-per-task=4 --mem=36G --time=24:00:00 \
  scripts/run-cap-appimage-slurm.sh \
  --appimage /shared/apps/CAP-<version>-x86_64.AppImage \
  --assembly /shared/data/genome.fasta \
  --outdir /shared/results/genome-cap
```

Partition, account, reservation, mail, and node-exclusion settings are not
hard-coded. Supply them as `sbatch` options or adapt the header for the site.
Optional inputs use `--te-gff`, `--gene-gff`, `--metadata`, `--templates`, and
`--trash2`. Repeat limits use `--max-rep-size` and `--min-rep-size`.

If compute nodes do not support FUSE, use the extraction procedure above on
shared storage and pass the extracted entry point instead:

```bash
sbatch scripts/run-cap-appimage-slurm.sh \
  --appimage /shared/apps/cap-appimage/squashfs-root/AppRun \
  --assembly /shared/data/genome.fasta \
  --outdir /shared/results/genome-cap
```

The wrapper creates a unique Nextflow work directory under
`${SLURM_TMPDIR:-${TMPDIR:-/tmp}}`. It retains the directory after both success
and failure and prints its exact path plus a cleanup command. Final output must
use persistent storage; node-local work may be lost when an allocation ends.
Keep it for diagnosis or Nextflow resume, and remove it manually when no longer
needed.

Run its scheduler-free behavioral tests with:

```bash
make test-appimage-slurm-wrapper
```

These tests use a fake AppImage to validate arguments, CPU propagation,
optional inputs, paths containing spaces, unique work-directory creation, and
retention without calling `sbatch`.
They cannot check your scheduler, compute-node FUSE support, shared filesystems,
or site policy. Ask your cluster administrators about those requirements before
submitting CAP.

## Offline use

The AppImage contains CAP and all required runtime dependencies. The pipeline
can run without downloading packages. Nextflow may print a harmless version
check warning when network access is unavailable.

For offline systems, download the AppImage and checksum on a connected machine,
verify them, and transfer both files. If the destination does not provide FUSE,
use the extraction procedure above.

## Published v1.0.3 artifact

The immutable Linux x86-64 AppImage release is:

```text
CAP-v1.0.3-x86_64.AppImage
Size: 636,085,440 bytes
SHA-256: c4b586e6dc1d73fd80cbdb6d7a141051164b8aafb72e9b8327328fcb477f12fc
Release: https://github.com/MellifluousJDG/CAP/releases/tag/v1.0.3
```

## Build from source

Building is supported on Linux x86-64. The build requires a recursive CAP Git
checkout, network access, and common host tools including Bash, Git, `curl`,
`readelf`, and `sha256sum`. Conda and a C++ compiler do not need to be installed:
the build downloads a checksum-verified Miniconda bootstrap locally.

```bash
version=v<version>
test -z "$(git status --porcelain --ignore-submodules=none)"
test "$(git describe --tags --exact-match)" = "$version"
APPIMAGE_VERSION="$version" make appimage
```

The build downloads the exact runtime packages named by
`conda-runtime-linux-64.lock`, compiles CTW/BCT with the pinned compiler and
sysroot from `conda-linux-64.lock`, and writes
`dist/CAP-${version}-x86_64.AppImage` and its checksum. Both temporary Conda
environments use only conda-forge and bioconda; the build does not consult
Anaconda's default channels. For non-release builds, omitting
`APPIMAGE_VERSION` retains the short-CAP-revision filename.

Test the release artifact with:

```bash
version=v<version>
make test-appimage APPIMAGE="dist/CAP-${version}-x86_64.AppImage"
```

This runs the complete one-core CAP regression in both normal FUSE-mounted
execution and from a separately extracted directory. Each run must produce 16
outputs, all 15 deterministic non-PNG files must match their expected SHA-256
values, and the PNG must be nonempty. It also checks key bundled runtimes and
confirms that excluded Windows dependency bundles are absent. Normal execution
therefore requires FUSE 2 support, including `libfuse.so.2`, on the test host.

## Publish a release

Publication is manual. Build and test from the exact annotated Git tag as shown
above. If the GitHub Release does not exist yet, create it and upload the
already-tested files without rebuilding or renaming them:

```bash
version=v<version>
gh release create "$version" \
  "dist/CAP-${version}-x86_64.AppImage" \
  "dist/CAP-${version}-x86_64.AppImage.sha256" \
  --repo MellifluousJDG/CAP \
  --verify-tag \
  --title "$version" \
  --notes "CAP $version"
```

If the GitHub Release already exists but does not yet contain these assets,
upload them without replacing any existing asset:

```bash
version=v<version>
gh release upload "$version" \
  "dist/CAP-${version}-x86_64.AppImage" \
  "dist/CAP-${version}-x86_64.AppImage.sha256" \
  --repo MellifluousJDG/CAP
```

GitHub CLI authentication must be configured interactively before either
command; do not place tokens on the command line or in shell history. Published
release tags and assets are immutable. Do not replace them.

Download both public assets into a new directory and verify the published
checksum:

```bash
version=v<version>
consumer_dir="../consumer-${version}"
test ! -e "$consumer_dir"
mkdir "$consumer_dir"
gh release download "$version" \
  --repo MellifluousJDG/CAP \
  --dir "$consumer_dir"
(
  cd "$consumer_dir"
  sha256sum --check "CAP-${version}-x86_64.AppImage.sha256"
)
```

Remove generated AppImage build directories and artifacts with:

```bash
make clean-appimage
```

## Troubleshooting

- `Permission denied`: run `chmod +x` on the AppImage. If the filesystem is
  mounted with `noexec`, move it to a user-writable executable filesystem or use
  extraction mode.
- FUSE or mount errors: use the extraction procedure above.
- Assembly not found: use an absolute path and check that the file is readable.
- Output or work directory is not writable: create it in a location owned by
  your user and pass its absolute path.
- Interrupted or failed run: inspect `.nextflow.log` in the launch directory
  and retain the work directory for diagnosis.
- Disk usage grows during a run: Nextflow intermediates are stored under the
  path supplied with `-work-dir`; remove them only after a successful run.
- Cluster execution: ask administrators whether AppImage/FUSE is supported on
  compute nodes. Extraction mode may work where FUSE is unavailable. Use
  persistent storage for final results and node-local temporary storage for the
  Nextflow work directory when available.
