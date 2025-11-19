#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ────────────────────── HELP MESSAGE ──────────────────────
if (!params.assembly) {
    def help = """
╔══════════════════════════════════════════════════════════════════════════════╗
║                                 CAP - Centromere Analysis Pipeline            ║
╚══════════════════════════════════════════════════════════════════════════════╝

Required:
  --assembly   <path>   Genome assembly in FASTA format

Optional:
  --te_gff     <path>   EDTA TE annotation (GFF3)
  --gene_gff   <path>   Helixer gene annotation (GFF)
  --templates  <path>   Custom repeat templates (FASTA)
  --cores      <int>    Number of CPU cores [default: ${params.cores}]
  --outdir     <path>   Output directory [default: ${params.outdir}]

────────────────────────────────────────────────────────────────────────────────
                             INSTALLATION & RUN OPTIONS
────────────────────────────────────────────────────────────────────────────────

1. DOCKER (Recommended – zero setup)
   docker run --rm -v \$(pwd)/data:/data yourname/cap-pipeline:latest \\
     nextflow run main.nf -profile docker --assembly /data/genome.fasta

   → No local installation needed
   → Works on HPC (Singularity), macOS, Windows (WSL)

2. CONDA (Local Linux/macOS/WSL)
   conda env create -f environment.yml
   conda activate cap-pipeline
   nextflow run main.nf -profile conda --assembly data/genome.fasta

   → Full dependency isolation
   → Great for development

3. R + renv (Lightweight, R-only)
   R -e 'renv::restore()'   # installs exact R packages
   nextflow run main.nf -profile renv --assembly data/genome.fasta

   → Only installs R packages
   → Requires R ≥ 4.3 and Nextflow installed

────────────────────────────────────────────────────────────────────────────────
                                   EXAMPLES
────────────────────────────────────────────────────────────────────────────────

# Docker (HPC or local)
docker run --rm -v \$(pwd)/test:/data yourname/cap-pipeline:latest \\
  nextflow run main.nf -profile docker \\
    --assembly /data/sample.fasta \\
    --te_gff /data/EDTA.gff3 \\
    --outdir /data/results

# Conda
nextflow run main.nf -profile conda \\
  --assembly data/genome.fasta \\
  --gene_gff data/helixer.gff

# renv
nextflow run main.nf -profile renv \\
  --assembly data/genome.fasta

────────────────────────────────────────────────────────────────────────────────
                                    OUTPUTS
────────────────────────────────────────────────────────────────────────────────
- *_CAP_plot.png
- *_CAP_dotplot.png
- *_CAP_repeat_families.csv
- *_centromeric_scores.csv
- *_predictions.csv
- *_metadata.csv

See README.md for full details: https://github.com/vlothec/CAP

"""
    println help
    System.exit(0)
}

// ────────────────────── PARAMETERS ──────────────────────
params.assembly      = null
params.templates     = null
params.te_gff        = null
params.gene_gff      = null
params.cores         = 4
params.max_rep_size  = 100000
params.min_rep_size  = 50
params.outdir        = "./results"

// ────────────────────── INPUT CHANNELS ──────────────────────
assembly_ch = Channel.fromPath(params.assembly, checkIfExists: true)

templates_ch = params.templates
    ? Channel.fromPath(params.templates)
    : Channel.empty()

te_gff_ch = params.te_gff
    ? Channel.fromPath(params.te_gff)
    : Channel.empty()

gene_gff_ch = params.gene_gff
    ? Channel.fromPath(params.gene_gff)
    : Channel.empty()

// ────────────────────── TRASH2 ──────────────────────
process TRASH2 {
    tag "TRASH2 on ${assembly}"
    cpus params.cores

    input:
    path assembly
    path templates

    output:
    tuple path("${assembly.baseName}_repeats_with_seq.csv"),
          path("${assembly.baseName}_arrays.csv") into trash_out_ch

    script:
    def t = templates.name != 'NO_FILE' ? "-t ${templates}" : ''
    """
    TRASH.R -f ${assembly} -o . \
        --cores_no ${task.cpus} \
        --max_rep_size ${params.max_rep_size} \
        --min_rep_size ${params.min_rep_size} \
        ${t}
    """
}

// ────────────────────── GET METADATA ──────────────────────
process GET_METADATA {
    input:
    path assembly from assembly_ch

    output:
    path "${assembly.baseName}_metadata.csv" into metadata_ch

    script:
    """
    Rscript ${baseDir}/bin/get_metadata.R ${assembly} ${assembly.baseName}_metadata.csv
    """
}

// ────────────────────── FILTER TRASH ──────────────────────
process FILTER_TRASH {
    input:
    tuple path(repeats), path(arrays) from trash_out_ch

    output:
    tuple path("${repeats.baseName}_filtered.csv"),
          path("${arrays.baseName}_filtered.csv") into filtered_trash_ch

    script:
    """
    Rscript ${baseDir}/bin/filter_trash.R ${repeats} ${arrays} \
        ${repeats.baseName}_filtered.csv ${arrays.baseName}_filtered.csv
    """
}

// ────────────────────── MERGE CLASSES ──────────────────────
process MERGE_CLASSES {
    input:
    tuple path(repeats_f), path(arrays_f) from filtered_trash_ch

    output:
    tuple path("${repeats_f.baseName}_reclassed.csv"),
          path("${arrays_f.baseName}_reclassed.csv"),
          path("${repeats_f.baseName}_genome_classes.csv") into merged_ch

    script:
    """
    Rscript ${baseDir}/bin/merge_classes.R ${repeats_f} ${arrays_f} \
        ${repeats_f.baseName}_reclassed.csv \
        ${arrays_f.baseName}_reclassed.csv \
        ${repeats_f.baseName}_genome_classes.csv
    """
}

// ────────────────────── OPTIONAL: PARSE / FILTER TEs ──────────────────────
process PARSE_TES {
    input:
    path te_gff from te_gff_ch

    output:
    path "${te_gff.baseName}_TEs_parsed.csv" into parsed_tes_ch

    when: te_gff.name != 'NO_FILE'

    script:
    """
    Rscript ${baseDir}/bin/parse_TEs.R ${te_gff} ${te_gff.baseName}_TEs_parsed.csv
    """
}

process FILTER_TES {
    input:
    path parsed from parsed_tes_ch

    output:
    path "${parsed.baseName}_filtered.csv" into filtered_tes_ch

    script:
    """
    Rscript ${baseDir}/bin/filter_TEs.R ${parsed} ${parsed.baseName}_filtered.csv
    """
}

// ────────────────────── OPTIONAL: PARSE / FILTER GENES ──────────────────────
process PARSE_GENES {
    input:
    path gene_gff from gene_gff_ch

    output:
    path "${gene_gff.baseName}_genes_parsed.csv" into parsed_genes_ch

    when: gene_gff.name != 'NO_FILE'

    script:
    """
    Rscript ${baseDir}/bin/parse_genes.R ${gene_gff} ${gene_gff.baseName}_genes_parsed.csv
    """
}

process FILTER_GENES {
    input:
    path parsed from parsed_genes_ch

    output:
    path "${parsed.baseName}_filtered.csv" into filtered_genes_ch

    script:
    """
    Rscript ${baseDir}/bin/filter_genes.R ${parsed} ${parsed.baseName}_filtered.csv
    """
}

// ────────────────────── SCORE CENTROMERIC ──────────────────────
process SCORE_CENTROMERIC {
    input:
    tuple path(repeats_r), path(arrays_r), path(genome_classes) from merged_ch
    path metadata from metadata_ch
    path te_f from filtered_tes_ch.mix(Channel.empty()).first()
    path gene_f from filtered_genes_ch.mix(Channel.empty()).first()

    output:
    path "${assembly.baseName}_centromeric_scores.csv" into scores_ch

    script:
    def te_opt  = te_f.name != 'NO_FILE' ? te_f : ''
    def gene_opt = gene_f.name != 'NO_FILE' ? gene_f : ''
    """
    Rscript ${baseDir}/bin/score_centromeric_classes.R \
        ${repeats_r} ${arrays_r} ${genome_classes} \
        ${metadata} \
        ${te_opt} ${gene_opt} \
        ${assembly.baseName}_centromeric_scores.csv
    """
}

// ────────────────────── PREDICT CENTROMERIC ──────────────────────
process PREDICT_CENTROMERIC {
    input:
    path scores from scores_ch

    output:
    path "${scores.baseName}_predictions.csv" into predictions_ch

    script:
    """
    python3 ${baseDir}/bin/predict_centromeric.py \
        ${scores} /model/centromeric_model_v2.pkl \
        ${scores.baseName}_predictions.csv
    """
}

// ────────────────────── GC & CTW (parallel) ──────────────────────
process GC {
    input:
    path assembly from assembly_ch

    output:
    path "${assembly.baseName}_GC.csv" into gc_ch

    script:
    """
    Rscript ${baseDir}/bin/GC.R ${assembly} ${assembly.baseName}_GC.csv
    """
}

process CTW {
    input:
    path assembly from assembly_ch

    output:
    path "${assembly.baseName}_CTW.csv" into ctw_ch

    script:
    """
    Rscript ${baseDir}/bin/CTW.R ${assembly} ${assembly.baseName}_CTW.csv
    """
}

// ────────────────────── FINAL CAP ──────────────────────
process CAP {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path predictions      from predictions_ch
    path repeats_r        from merged_ch.map { it[0] }
    path arrays_r         from merged_ch.map { it[1] }
    path genome_classes   from merged_ch.map { it[2] }
    path metadata         from metadata_ch
    path te_f             from filtered_tes_ch.mix(Channel.empty()).first()
    path gene_f           from filtered_genes_ch.mix(Channel.empty()).first()

    output:
    path "*_CAP_plot.png"
    path "*_CAP_dotplot.png"
    path "*_CAP_repeat_families.csv"
    path "*_CAP_model.txt"
    path "*_TEs_filtered.csv"   optional true
    path "*_genes_filtered.csv" optional true

    script:
    def te_opt  = te_f.name != 'NO_FILE' ? te_f : ''
    def gene_opt = gene_f.name != 'NO_FILE' ? gene_f : ''
    """
    Rscript ${baseDir}/bin/CAP.R \
        ${predictions} ${repeats_r} ${arrays_r} ${genome_classes} \
        ${metadata} \
        ${te_opt} ${gene_opt} \
        ${assembly.baseName}_CAP_plot.png \
        ${assembly.baseName}_CAP_repeat_families.csv \
        ${assembly.baseName}_CAP_dotplot.png \
        ${assembly.baseName}_CAP_model.txt
    """
}

// ────────────────────── WORKFLOW ──────────────────────
workflow {
    // ---- TRASH2 chain ----
    TRASH2(assembly_ch, templates_ch)
    FILTER_TRASH()
    MERGE_CLASSES()

    // ---- Optional TEs ----
    if (params.te_gff) {
        PARSE_TES()
        FILTER_TES()
    }

    // ---- Optional Genes ----
    if (params.gene_gff) {
        PARSE_GENES()
        FILTER_GENES()
    }

    // ---- Scoring + prediction ----
    SCORE_CENTROMERIC()
    PREDICT_CENTROMERIC()

    // ---- GC / CTW in parallel ----
    GC()
    CTW()

    // ---- Final visualisation ----
    CAP()
}