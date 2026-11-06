# CAP: Centromere Analysis Pipeline

**CAP** is a Nextflow pipeline for centromere prediction and repeat analysis in genomic assemblies. It integrates **TRASH_2** for tandem-repeat detection, processes optional transposon and gene annotations, and produces plots and tables.

## Published v1.0.3 artifacts

The current immutable Linux release artifacts are:

```text
AppImage: CAP-v1.0.3-x86_64.AppImage
SHA-256: c4b586e6dc1d73fd80cbdb6d7a141051164b8aafb72e9b8327328fcb477f12fc

OCI: ghcr.io/mellifluousjdg/cap:1.0.3
Digest: sha256:2563036a264a3d7b72bba57215229df599af15a36ae0c4d7155b70a3127b6f95

Source: ghcr.io/mellifluousjdg/cap-source:1.0.3
Digest: sha256:33daf3e3b1b90011201c528d93b3bd86052258cf9a410bb1d1376c4a2a0699e0
```

See `docs/OCI_RELEASE.md` for the source archive size and SHA-256 and for the
manual publication and verification procedure.

## AppImage quick start

For Linux x86-64, CAP v1.0.3 is available from the GitHub Release page as
`CAP-v1.0.3-x86_64.AppImage` (636,085,440 bytes, SHA-256
`c4b586e6dc1d73fd80cbdb6d7a141051164b8aafb72e9b8327328fcb477f12fc`).
Download the AppImage and its `.sha256` file. The AppImage contains CAP and all
runtime dependencies, so Conda and a container runtime are not required.

```bash
sha256sum --check CAP-<version>-x86_64.AppImage.sha256
chmod +x CAP-<version>-x86_64.AppImage
mkdir -p /path/to/cap-results /path/to/cap-work

./CAP-<version>-x86_64.AppImage \
  --assembly /path/to/genome.fasta \
  --cores 1 \
  --outdir /path/to/cap-results \
  -work-dir /path/to/cap-work
```

Replace `<version>` with the version shown in the release. If the system does
not support FUSE, extract the AppImage with `--appimage-extract` and run
`./squashfs-root/AppRun` instead.

See `docs/APPIMAGE.md` for download verification, normal and extracted
execution, optional inputs, offline use, building from source, testing, cleanup,
and troubleshooting.

## Conda quick start

The tested packaging path is Conda on Linux x86-64. Conda must already be installed; Nextflow, Java, R, Python, Git, and scientific dependencies are installed in one CAP environment.

```bash
git clone --recurse-submodules \
  https://github.com/MellifluousJDG/CAP.git
cd CAP

# Create the default cap-pipeline environment and compile CTW/BCT.
make install

# Run the bundled one-core regression test and verify all 16 checksums.
make test-conda

# Run CAP on a genome assembly in FASTA format.
make run ASSEMBLY=/path/to/genome.fasta
```

`ASSEMBLY` is a Make variable containing the path to the input genome assembly. The `run` target passes it to Nextflow as `--assembly`.

## Conda installation

### Reproducible Linux installation

On Linux x86-64, `setup_conda.sh` creates the environment from the committed `conda-linux-64.lock`. The lock contains the exact package URLs used by the tested environment.

On other platforms, setup falls back to solving `environment.yml`. This is a
best-effort installation path without the exact-package reproducibility or test
coverage of Linux x86-64.

A recursive clone is preferred:

```bash
git clone --recurse-submodules \
  https://github.com/MellifluousJDG/CAP.git
cd CAP
make install
```

If CAP was cloned without `--recurse-submodules`, `make install` initializes the TRASH_2 submodule using Git from the Conda environment.

The setup script does not activate the environment in the calling shell. To work in it interactively, run:

```bash
conda activate cap-pipeline
```

The Make targets use `conda run`, so activation is not required for `make run` or `make test-conda`.

### Custom environment name

The default environment name is `cap-pipeline`. Use `ENV_NAME` consistently with Make:

```bash
make install ENV_NAME=my-cap
make test-conda ENV_NAME=my-cap
make run ENV_NAME=my-cap ASSEMBLY=/path/to/genome.fasta
make remove-env ENV_NAME=my-cap
```

When invoking the setup script directly, use `CAP_ENV_NAME`:

```bash
CAP_ENV_NAME=my-cap bash setup_conda.sh
```

### Make targets

The Makefile is a curated convenience interface for common installation,
execution, testing, packaging, and maintenance tasks; it is not an exhaustive
list of every supported command. Specialized release, compliance, and cluster
operations are invoked directly as documented in the relevant guides under
`docs/`.

- `make install` creates or reuses the environment, initializes TRASH_2 when necessary, and compiles CTW/BCT.
- `make compile` recompiles CTW/BCT in the selected environment.
- `make run ASSEMBLY=...` runs CAP on the supplied FASTA assembly.
- `make run-trash2 ASSEMBLY=... OUTDIR=...` runs only the pinned TRASH_2 module.
- `make lint` checks every tracked shell and Python script for syntax errors.
- `make test-ctw-build` runs isolated fixtures against the shared CTW build and execution validator; it does not require a Conda environment or modify `bin/ctw-calc`.
- `make test-conda` performs the bundled deterministic one-core regression test.
- `make test-trash2-standalone` validates the standalone launcher without running a full analysis.
- `make remove-env` removes the selected Conda environment.
- `make clean-build` removes the generated `bin/ctw-calc` binary without removing the environment.

### Removing the environment

```bash
make remove-env
```

This removes the Conda environment but not Nextflow work or result directories.

## Usage

Run through Make:

```bash
make run ASSEMBLY=/path/to/genome.fasta
```

For options beyond `ASSEMBLY` and `ENV_NAME`, invoke Nextflow directly as shown below.

Or activate the environment and invoke Nextflow directly:

```bash
conda activate cap-pipeline
nextflow run . --assembly /path/to/genome.fasta
```

### Standalone TRASH_2

After `make install`, you can run only CAP's pinned TRASH_2 module with the same locked
Conda environment:

```bash
make run-trash2 \
  ASSEMBLY=/path/to/genome.fasta \
  OUTDIR=/path/to/trash2-results
```

The output directory is created when its parent exists. Common options are:

```bash
make run-trash2 \
  ASSEMBLY=/path/to/genome.fasta \
  OUTDIR=/path/to/trash2-results \
  ENV_NAME=cap-pipeline \
  CORES=1 \
  MAX_REP_SIZE=1000 \
  MIN_REP_SIZE=7 \
  TEMPLATES=/path/to/repeat-templates.fasta
```

`TEMPLATES` is optional. `ENV_NAME` must identify an existing environment made
by `make install`; this command never creates or updates an environment. It
checks the locked R packages, MAFFT, and HMMER before starting and fails if a
dependency is absent instead of installing packages at runtime.

For an assembly named `genome.fasta`, principal outputs include:

- `genome.fasta_repeats_with_seq.csv`;
- `genome.fasta_repeats.csv` and `genome.fasta_repeats.gff`;
- `genome.fasta_arrays.csv` and `genome.fasta_arrays.gff`;
- `genome.fasta_run_time.csv`.

Run `scripts/run-trash2.sh --help` for the equivalent direct launcher options.

### Required argument

- `--assembly`: path to the genome assembly in FASTA format.

A genome assembly contains one or more reconstructed genomic DNA sequences, such as chromosomes, scaffolds, or contigs.

### Optional arguments

- `--te_gff`: EDTA transposable-element annotation in GFF3 format; default `null`.
- `--gene_gff`: Helixer gene annotation in GFF format; default `null`.
- `--metadata`: existing chromosome metadata CSV; default `null`.
- `--templates`: repeat-template FASTA passed to TRASH_2; default `null`.
- `--trash2`: existing TRASH_2 results directory, which skips the TRASH_2 run; default `null`.
- `--cores`: number of CPU cores; default `1`.
- `--outdir`: output directory; default `./results`.
- `--max_rep_size`: maximum tandem-repeat size; default `1000`.
- `--min_rep_size`: minimum tandem-repeat size; default `7`.

Example with annotations:

```bash
conda activate cap-pipeline
nextflow run . \
  --assembly data/arabidopsis.fasta \
  --te_gff data/EDTA.TEanno.split.gff3 \
  --gene_gff data/helixer.gff \
  --cores 1 \
  --outdir results_arabidopsis
```

## Reproducibility and smoke testing

Run:

```bash
make test-conda
```

The test:

1. uses the bundled small FASTA assembly;
2. requires the TRASH_2 submodule and compiled `bin/ctw-calc`;
3. runs the complete ten-process workflow with one core and fresh result/work directories;
4. requires 16 nonempty output files, including CSV, PNG, and RDS output;
5. verifies every result against `test/expected-results.sha256`.

The exact checksum guarantee currently applies to the tested Linux x86-64 environment, committed model and dependencies, bundled input, published TRASH_2 revision, and one-core execution. A fixed sampling seed is used, but deterministic multiprocessing has not yet been established. Do not interpret the one-core regression test as a guarantee that every multicore run is byte-identical.

`test/ath_Chr1_extraction_trc.fasta` is the committed regression input, and
`test/expected-results.sha256` records the expected output hashes. Generated
result and Nextflow work files remain uncommitted.

Override disposable test locations if needed:

```bash
make test-conda \
  RESULTS_DIR=/tmp/cap-results \
  WORK_DIR=/tmp/cap-work
```

## Updating dependencies and the Linux lock

`environment.yml` is the human-maintained dependency specification. `conda-linux-64.lock` is the exact tested Linux x86-64 artifact and should not be edited casually.

When dependencies must change:

1. edit `environment.yml` intentionally;
2. solve a fresh Linux x86-64 environment using only its declared channels;
3. test imports and model loading, including scikit-learn and XGBoost compatibility;
4. compile CTW/BCT;
5. run `make test-conda` and investigate every output difference;
6. export exact package URLs from the validated environment:

   ```bash
   conda list -n cap-pipeline --explicit > conda-linux-64.lock
   ```

7. from a fresh recursive clone containing the candidate lock, create a
   disposable environment from that lock and run the complete smoke test:

   ```bash
   conda create -n cap-lock-test --file conda-linux-64.lock
   make compile ENV_NAME=cap-lock-test
   make test-conda ENV_NAME=cap-lock-test
   ```

8. update `test/expected-results.sha256` only when output changes are understood and approved;
9. commit dependency specifications, lock changes, and intentional output-baseline changes separately where practical.

The lock must not contain local filesystem paths or private package URLs. Review the generated file before committing it.

## Outputs

CAP writes key files such as:

- `*_CAP_plot*.png`: centromere, repeat-density, and annotation visualization;
- `*_CAP_model.txt`: text summary of identified centromeres;
- `*_CAP_repeat_families.csv`: repeat-family statistics;
- `*_centromeric_scores.csv`: calculated centromeric features and scores;
- `*_centromeric_scores_predictions.csv`: model predictions;
- `*_CAP_Rdata.rds`: serialized R data used by the final report.

Additional intermediate CSV files contain detected, filtered, and reclassified repeats.

## Pipeline overview

CAP performs six main stages:

1. **Repeat identification:** run TRASH_2 to identify tandem repeats, or reuse validated precomputed TRASH_2 results.
2. **Filtering and merging:** filter detected repeats and merge related repeat classes.
3. **Feature extraction:** calculate sequence features and process optional TE and gene annotations.
4. **Scoring:** combine repeat, sequence, and annotation information to score candidate centromeric regions.
5. **Prediction:** apply the bundled machine-learning model to the calculated scores.
6. **Visualization and reporting:** produce plots, summaries, tables, and reusable R data.

In the diagram below, `A ──> B` means that B waits for and uses output from A.
Branches on separate lines may run at the same time, while joined arrows mean that the next stage waits for all required inputs.

```text
Genome assembly
    |
    +──> Identify, filter, and classify repeats ──+
    |                                             |
    +──> Prepare metadata ────────────────────────+──> Score candidates
    |                                             |         |
    +──> Process optional annotations ────────────+         v
    |                                                   Predict
    +──> Calculate GC and CTW ──────────────────────────+  |
                                                        |  |
                                                        v  v
                                                  Plots and reports
```

Nextflow may run independent stages at the same time and waits automatically when one stage needs another stage's output.
See `docs/PIPELINE_FLOW.md` for the complete technical process graph, filename rules, output inventory, precomputed TRASH_2 behavior, and troubleshooting guidance.

## Repository structure

```text
CAP/
├── main.nf                         # Nextflow workflow
├── nextflow.config                 # Parameters and execution profiles
├── Makefile                        # Local commands, tests, and packaging entry points
├── setup_conda.sh                  # Locked Conda setup and CTW/BCT compilation
├── environment.yml                 # Human-maintained development dependencies
├── conda-linux-64.lock             # Exact tested development packages
├── environment-runtime.yml         # Human-maintained OCI runtime dependencies
├── conda-runtime-linux-64.lock     # Exact tested OCI runtime packages
├── Dockerfile                      # Reproducible OCI image definition
├── scripts/
│   ├── lint.sh                     # Shell and Python syntax checks
│   ├── build-ctw.sh                # Shared CTW build and execution validation
│   ├── test-ctw-build.sh           # Isolated CTW build-validator fixtures
│   ├── cap-container-entrypoint.sh # OCI runtime entry point
│   ├── run-cap-slurm.sh            # Conda single-allocation SLURM launcher
│   ├── run-cap-appimage-slurm.sh   # AppImage single-allocation SLURM launcher
│   ├── run-trash2.sh               # Locked standalone TRASH_2 launcher
│   ├── test-conda.sh               # One-core CAP regression test
│   ├── test-container.sh           # OCI lock and compliance checks
│   ├── test-container-image.sh     # Built OCI image integration test
│   ├── test-slurm-wrapper.sh       # Scheduler-free Conda SLURM launcher tests
│   ├── test-appimage-slurm-wrapper.sh # Scheduler-free AppImage SLURM tests
│   ├── test-trash2-standalone.sh   # Standalone TRASH_2 launcher tests
│   ├── appimage/
│   │   ├── AppRun                  # AppImage entry point
│   │   ├── Rscript                 # AppImage Rscript path wrapper
│   │   ├── build-appimage.sh       # AppImage build script
│   │   └── test-built-appimage.sh  # Built AppImage regression test
│   └── compliance/
│       ├── export-debian-compliance-inputs.sh # Export Debian package metadata
│       ├── finalize-compliance-checksums.sh   # Write bundle checksums
│       ├── generate-oci-compliance-bundle.sh # Build OCI compliance bundle
│       ├── generate-oci-inventory.sh         # Export OCI package inventory
│       ├── test-collect-conda-notices.sh     # Conda notice collector fixtures
│       ├── test-collect-conda-sources.sh     # Conda source collector fixtures
│       ├── test-collect-debian-sources.py    # Debian source collector fixtures
│       ├── test-conda-source-notices.sh      # Source-notice review fixtures
│       ├── test-finalize-compliance-checksums.sh # Bundle checksum fixtures
│       ├── test-merge-conda-source-corpora.py # Source-corpus merge fixtures
│       ├── test-sbom.sh                      # CycloneDX SBOM fixtures
│       ├── test-source-review-policy.py      # Source-review policy fixtures
│       └── test-source-obligations.sh        # License-obligation fixtures
├── docs/                           # Usage, architecture, and release guides
├── bin/                            # Pipeline scripts and CTW source
├── modules/TRASH_2/                # Pinned TRASH_2 Git submodule
├── model/                          # Pre-trained prediction models
└── test/
    ├── ath_Chr1_extraction_trc.fasta # Regression-test assembly input
    └── expected-results.sha256       # Deterministic output checksums
```

Run `make lint` before committing changes to shell or Python scripts. It checks
all tracked scripts, including extensionless AppImage entry points, without
writing bytecode into the repository.

## Single-node SLURM execution

CAP includes `scripts/run-cap-slurm.sh` for clusters where one allocation should run the entire workflow. SLURM starts this script once, and Nextflow uses its local executor inside the allocated node. Nextflow does not submit additional scheduler jobs, so dependent stages do not return to the queue.

Install CAP and its shared Conda environment on storage visible to the compute nodes, then submit:

```bash
sbatch scripts/run-cap-slurm.sh \
  --assembly /shared/data/genome.fasta \
  --outdir /shared/results/genome-cap
```

The script uses `SLURM_CPUS_PER_TASK` for `--cores`. Its generic defaults request one CPU and 36 GB of memory. Override resources with `sbatch` options when submitting, for example:

```bash
sbatch --cpus-per-task=4 --mem=36G --time=24:00:00 \
  scripts/run-cap-slurm.sh \
  --assembly /shared/data/genome.fasta \
  --outdir /shared/results/genome-cap
```

Partition, account, excluded nodes, and site-specific paths are deliberately not hard-coded. Add them as `sbatch` options or adapt the `#SBATCH` header for the target cluster.

The default Conda environment is `cap-pipeline`. Select another named environment or Conda executable with:

```bash
sbatch scripts/run-cap-slurm.sh \
  --assembly /shared/data/genome.fasta \
  --outdir /shared/results/genome-cap \
  --env-name my-cap \
  --conda /shared/tools/miniconda3/bin/conda
```

Optional annotations use `--te-gff`, `--gene-gff`, and `--metadata`; repeat
templates use `--templates`. Repeat limits use `--max-rep-size` and
`--min-rep-size`. To reuse existing TRASH_2 output explicitly:

```bash
sbatch scripts/run-cap-slurm.sh \
  --assembly /shared/data/genome.fasta \
  --outdir /shared/results/genome-cap \
  --trash2 /shared/results/precomputed-trash2
```

For an assembly named `genome.fasta`, that directory must contain both `genome.fasta_repeats_with_seq.csv` and `genome.fasta_arrays.csv`. A supplied but invalid `--trash2` directory causes an error instead of silently rerunning TRASH_2. Omit `--trash2` to run TRASH_2 normally.

The wrapper creates a unique Nextflow work directory under
`${SLURM_TMPDIR:-${TMPDIR:-/tmp}}`, retains it after both successful and failed
runs, and prints its exact path plus a cleanup command. Final output must use
persistent/shared storage. Keep the work directory for diagnosis or Nextflow
resume; remove it manually when it is no longer needed.

The script refuses to run directly on a login node. Its explicit local-validation mode is:

```bash
CAP_SLURM_LOCAL_TEST=1 bash scripts/run-cap-slurm.sh \
  --assembly test/ath_Chr1_extraction_trc.fasta \
  --outdir results_slurm_local
```

Run the scheduler-free launcher tests with:

```bash
make test-slurm-wrapper
```

These tests validate argument handling, CPU propagation, optional TRASH_2
behavior, unique work-directory creation, and unconditional retention without
calling `sbatch`.

For AppImage-based clusters, use `scripts/run-cap-appimage-slurm.sh` instead.
See `docs/APPIMAGE.md` for submission, extracted-AppImage use, and the limits of
scheduler-free validation.

## OCI container quick start

The public Linux x86-64 image contains CAP and all runtime dependencies. With
rootless Podman installed, put a FASTA file in a new working directory:

```bash
assembly=/path/to/genome.fasta
mkdir -p cap-run/data cap-run/results cap-run/work
cp "$assembly" cap-run/data/
assembly_name=$(basename "$assembly")
cd cap-run

image=ghcr.io/mellifluousjdg/cap:1.0.3
podman run --rm --userns=keep-id \
  -v "$PWD/data:/data:ro" \
  -v "$PWD/results:/results" \
  -v "$PWD/work:/work" \
  "$image" \
  --assembly "/data/$assembly_name" \
  --cores 1 \
  --outdir /results \
  -work-dir /work/nextflow
```

CAP writes final
CSV, RDS, model, and PNG files to `results/`; `work/nextflow` contains disposable
Nextflow intermediates. Successful rootless Podman runs leave both directories
owned by the invoking user. Remove `work/` after success to reclaim space.

For an immutable pull, replace the image assignment with:

```bash
image=ghcr.io/mellifluousjdg/cap@sha256:2563036a264a3d7b72bba57215229df599af15a36ae0c4d7155b70a3127b6f95
```

See `docs/OCI_IMAGE.md` for Docker, permissions, SELinux, troubleshooting,
cleanup, image details, and the published corresponding-source artifact; see
`docs/OCI_HPC.md` for single-allocation cluster templates.

## License

CAP's original code is distributed under the [MIT License](LICENSE). The
BCT-derived source in `bin/src/BCT` and its compiled `bin/ctw-calc` executable
use `GPL-2.0-or-later`, based on the `GPL (>= 2)` CRAN metadata for the related
BCT package published by the same author. See `bin/src/BCT/README.md` for the
license basis, provenance, build instructions, and corresponding-source
requirements.

## Contact

For questions or issues, open an issue on GitHub.
