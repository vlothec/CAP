# CAP pipeline flow and file naming

This document describes the data flow implemented in `main.nf`, the files
published to `--outdir`, and how optional inputs change execution.

## Inputs

CAP requires:

- `--assembly`: a genome assembly in FASTA format.

CAP optionally accepts:

- `--templates`: repeat templates passed to TRASH_2;
- `--te_gff`: EDTA transposable-element annotation;
- `--gene_gff`: Helixer gene annotation;
- `--metadata`: an existing chromosome metadata CSV used by `GET_METADATA`;
- `--trash2`: a directory containing precomputed TRASH_2 output;
- `--cores`: CPU count used by process resource configuration and passed to TRASH_2;
- `--max_rep_size`: maximum repeat size passed to TRASH_2;
- `--min_rep_size`: minimum repeat size passed to TRASH_2;
- `--outdir`: directory to which process outputs are copied.

Annotation sequence names and coordinates must correspond to the supplied assembly.

## End-to-end flow

In the diagram below, `A ──> B` means that B waits for and consumes output from
A. Separate branches may run concurrently. A stage with several incoming arrows
waits for every required input.

```text
assembly FASTA
    │
    ├──> CHECK_DEPS ──> TRASH2 ──────────────────────────┐
    │                   (or validated --trash2 files)    │
    │                                                    ▼
    │                                              FILTER_TRASH
    │                                                    │
    │                                                    └──> MERGE_CLASSES ──┐
    │                                                                         │
    ├──> GET_METADATA ────────────────────────────────────────────────────────┤
    │                                                                         │
    ├──> optional --te_gff ──> PARSE_TES ──> FILTER_TES ──────────────────────┼──> SCORE_CENTROMERIC
    │                                                                         │            │
    ├──> optional --gene_gff ──> PARSE_GENES ──> FILTER_GENES ────────────────┘            ▼
    │                                                                              PREDICT_CENTROMERIC
    ├──> GC ─────────────────────────────────────────────────────────────────────────────┐  │
    │                                                                                    │  │
    └──> CTW ────────────────────────────────────────────────────────────────────────────┤  │
                                                                                         ▼  ▼
                                                                                          CAP
                                                                                           │
                                                                                           ├──> plots
                                                                                           ├──> repeat-family table
                                                                                           ├──> model summary
                                                                                           └──> R data
```

When neither optional annotation is supplied, CAP skips that parsing branch and
continues without the corresponding annotation files. `GC` and `CTW` do not feed
`SCORE_CENTROMERIC`; they join the other results at the final `CAP` process.

Nextflow may execute independent branches concurrently. For example, metadata,
GC, and CTW do not need to wait for repeat filtering. Processes that consume
upstream output wait for that output automatically.

## Filename rules

Two assembly-derived forms appear in output names:

- **assembly name**: the complete input filename, including its final extension;
- **assembly stem**: the filename with its final extension removed.

For an input named `genome.fasta`:

```text
assembly name = genome.fasta
assembly stem = genome
```

TRASH_2 output retains the complete assembly name:

```text
genome.fasta_repeats_with_seq.csv
genome.fasta_arrays.csv
```

Most downstream assembly-derived output uses the stem:

```text
genome_metadata.csv
genome_GC.csv
genome_CTW.csv
genome_centromeric_scores.csv
genome_CAP_model.txt
```

Each downstream process appends its suffix to the stem of its immediate input.
Consequently, repeat filenames accumulate processing-stage suffixes such as
`_filtered` and `_reclassed`.

## Process and output reference

### `CHECK_DEPS`

Checks that `nhmmer` and `mafft` are available before repeat detection. It
publishes no file.

### `TRASH2`

Runs `modules/TRASH_2/src/TRASH.R` on the assembly. It accepts optional
templates and uses `--cores`, `--max_rep_size`, and `--min_rep_size`.

Published files:

```text
<assembly-name>_repeats_with_seq.csv
<assembly-name>_arrays.csv
```

The first file records repeat information and representative sequences. The
second records tandem-repeat arrays and their genomic coordinates.

### `FILTER_TRASH`

Filters the two TRASH_2 tables.

Published files:

```text
<assembly-name>_repeats_with_seq_filtered.csv
<assembly-name>_arrays_filtered.csv
```

### `MERGE_CLASSES`

Merges related repeat classes and produces reclassified repeat and array tables plus a genome-class summary.

Published files:

```text
<assembly-name>_repeats_with_seq_filtered_reclassed.csv
<assembly-name>_arrays_filtered_reclassed.csv
<assembly-name>_repeats_with_seq_filtered_genome_classes.csv
```

### `GET_METADATA`

Builds the chromosome metadata used downstream. When `--metadata` is supplied,
`get_metadata.R` uses that file while producing CAP's normalized metadata
output.

Published file:

```text
<assembly-stem>_metadata.csv
```

### `PARSE_TES` and `FILTER_TES`

These processes run only when `--te_gff` is supplied. They parse TE annotation
and filter it using array and chromosome context.

Published files:

```text
<TE-stem>_TEs_parsed.csv
<TE-stem>_TEs_parsed_filtered.csv
```

Without `--te_gff`, CAP skips these processes and publishes neither TE file.

### `PARSE_GENES` and `FILTER_GENES`

These processes run only when `--gene_gff` is supplied. They parse gene
annotation and filter it using array and chromosome context.

Published files:

```text
<gene-stem>_genes_parsed.csv
<gene-stem>_genes_parsed_filtered.csv
```

Without `--gene_gff`, CAP skips these processes and publishes neither gene file.

### `SCORE_CENTROMERIC`

Combines reclassified repeats, metadata, and optional TE/gene context into model features and centromeric scores.

Published file:

```text
<assembly-stem>_centromeric_scores.csv
```

### `PREDICT_CENTROMERIC`

Loads `model/centromeric_model_v2.pkl` and applies it to the score table.

Published file:

```text
<assembly-stem>_centromeric_scores_predictions.csv
```

### `GC`

Calculates genome sequence GC content for the final analysis.

Published file:

```text
<assembly-stem>_GC.csv
```

### `CTW`

Calculates the context-tree-weighting sequence metric used by the final analysis.

Published file:

```text
<assembly-stem>_CTW.csv
```

### `CAP`

Combines predictions, repeat classifications, metadata, optional annotations,
GC, CTW, and scores to produce the final visualization and reports.

Published files:

```text
<assembly-stem>_CAP_plot_*.png
<assembly-stem>_CAP_repeat_families.csv
<assembly-stem>_CAP_model.txt
<assembly-stem>_CAP_Rdata.rds
```

The plot filename records annotation availability. For example, the bundled
unannotated fixture produces a suffix containing `noedta_nogene`.

## Reusing precomputed TRASH_2 output

When `--trash2 DIR` is supplied, CAP skips the `TRASH2` process and requires these exact files in `DIR`:

```text
DIR/<assembly-name>_repeats_with_seq.csv
DIR/<assembly-name>_arrays.csv
```

For `genome.fasta`, the expected files are:

```text
DIR/genome.fasta_repeats_with_seq.csv
DIR/genome.fasta_arrays.csv
```

If either file is absent, CAP fails instead of silently running TRASH_2. The
supplied files feed `FILTER_TRASH` directly.

Because the `TRASH2` process is skipped, CAP does not copy those two source
files into the new `--outdir`; they remain in the supplied directory. The
downstream filtered, reclassified, scoring, prediction, and final files are
still published normally. Thus the bundled unannotated run publishes 16 files
when it runs TRASH_2 and 14 newly generated files when it reuses both TRASH_2
files.

## Published output versus Nextflow work

Each file-producing process uses `publishDir` in copy mode:

- task execution and staged files live under Nextflow's work directory;
- user-facing copies are written to `--outdir`;
- deleting a work directory removes task caches and diagnostics, not successfully published copies;
- incomplete or failed tasks may leave useful logs and temporary files only in their work directory.

The single-node SLURM wrapper creates a unique work directory under node-local
temporary storage while keeping `--outdir` on persistent storage. It retains the
work directory after every run and prints its exact path and a manual cleanup
command.

## Diagnosing an incomplete run

Nextflow's terminal output and `.nextflow.log` show the last process states. A
missing downstream file often means an upstream dependency failed rather than
that the final process independently omitted it.

Useful checks include:

1. find the first process marked failed in `.nextflow.log`;
2. inspect that task's `.command.sh`, `.command.out`, and `.command.err` in the reported work directory;
3. confirm the assembly and optional annotations use corresponding sequence names;
4. when using `--trash2`, confirm both exact assembly-name-prefixed files exist;
5. distinguish the persistent `--outdir` from Nextflow's task work directory.

The bundled deterministic fixture's complete expected filename set is recorded in `test/expected-results.sha256`.
