#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



// ────────────────────── PARAMETERS ──────────────────────
params.assembly = null
params.templates = null
params.te_gff = null
params.gene_gff = null
params.cores = 1
params.max_rep_size = 200
params.outdir = "./results"

// ────────────────── CHECK DEPENDANCIES ──────────────────
process CHECK_DEPS {
    tag "check nhmmer and mafft availability"
	output:
    val true
    script:
     """
    echo "Checking for nhmmer:"
    if ! command -v nhmmer &> /dev/null; then
        echo "❌ nhmmer not found in PATH"; exit 1;
    else
        nhmmer -h 2>&1 | head -n 3 || true
    fi

    echo "Checking for mafft:"
    if ! command -v mafft &> /dev/null; then
        echo "❌ mafft not found in PATH"; exit 1;
    else
        mafft --version 2>&1 || true
    fi

    echo "✅ All dependencies found."
    """
}



// ──────────────────────── TRASH2 ────────────────────────
process TRASH2 {
	publishDir "${params.outdir}", mode: 'copy'
    tag "TRASH2 on ${assembly.baseName}"
    cpus params.cores

    input:
    tuple (path(assembly), val(templates))
	val ready

	when:
	ready
	
    output:
    tuple (path("${assembly.name}_repeats_with_seq.csv"),
          path("${assembly.name}_arrays.csv"))
		  
    script:
	def t = params.templates ? "-t ${params.templates}" : ''
    """
    echo "Running TRASH2 with output dir: \$PWD"
    ${workflow.projectDir}/modules/TRASH_2/src/TRASH.R -f ${assembly} -o \$PWD \
        --cores_no ${task.cpus} \
        --max_rep_size ${params.max_rep_size} \
        ${t}
    """
}

// ────────────────────── GET METADATA ──────────────────────
process GET_METADATA {
	publishDir "${params.outdir}", mode: 'copy'
    input:
    path assembly
    output:
    path "${assembly.baseName}_metadata.csv"
    script:
    """
    Rscript ${workflow.projectDir}/bin/get_metadata.R ${assembly} ${assembly.baseName}_metadata.csv
    """
}

// ────────────────────── FILTER TRASH ──────────────────────
process FILTER_TRASH {
    input:
    tuple path(repeats), path(arrays)
    output:
    tuple path("${repeats.baseName}_filtered.csv"),
          path("${arrays.baseName}_filtered.csv")
    script:
    """
    Rscript ${workflow.projectDir}/bin/filter_trash.R ${repeats} ${arrays} \
        ${repeats.baseName}_filtered.csv ${arrays.baseName}_filtered.csv
    """
}

// ────────────────────── MERGE CLASSES ──────────────────────
process MERGE_CLASSES {
	publishDir "${params.outdir}", mode: 'copy'
    input:
    tuple path(repeats_f), path(arrays_f)
    output:
    tuple path("${repeats_f.baseName}_reclassed.csv"),
          path("${arrays_f.baseName}_reclassed.csv"),
          path("${repeats_f.baseName}_genome_classes.csv")
    script:
    """
    Rscript ${workflow.projectDir}/bin/merge_classes.R ${repeats_f} ${arrays_f} \
        ${repeats_f.baseName}_reclassed.csv \
        ${arrays_f.baseName}_reclassed.csv \
        ${repeats_f.baseName}_genome_classes.csv
    """
}

// ────────────────────── OPTIONAL: PARSE / FILTER TEs ──────────────────────
process PARSE_TES {
    input:
    path te_gff
    output:
    path "${te_gff.name}_TEs_parsed.csv"
	
    script:
    """
    Rscript ${workflow.projectDir}/bin/parse_TEs.R ${te_gff} ${te_gff.baseName}_TEs_parsed.csv
    """
}

process FILTER_TES {
    input:
    path parsed
    output:
    path "${parsed.baseName}_filtered.csv"
    script:
    """
    Rscript ${workflow.projectDir}/bin/filter_TEs.R ${parsed} ${parsed.baseName}_filtered.csv
    """
}

// ────────────────────── OPTIONAL: PARSE / FILTER GENES ──────────────────────
process PARSE_GENES {
    input:
    path gene_gff
    output:
    path "${gene_gff.name}_genes_parsed.csv"
 
    script:
    """
    Rscript ${workflow.projectDir}/bin/parse_genes.R ${gene_gff} ${gene_gff.baseName}_genes_parsed.csv
    """
}

process FILTER_GENES {
    input:
    path parsed
    output:
    path "${parsed.baseName}_filtered.csv"
    script:
    """
    Rscript ${workflow.projectDir}/bin/filter_genes.R ${parsed} ${parsed.baseName}_filtered.csv
    """
}

// ────────────────────── SCORE CENTROMERIC ──────────────────────
process SCORE_CENTROMERIC {
    publishDir "${params.outdir}", mode: 'copy'
	
    input:
    path assembly
    tuple path(repeats_r), path(arrays_r), path(genome_classes)
    path metadata
    val te_f
    val gene_f
	
    output:
    path "${assembly.baseName}_centromeric_scores.csv"

    script:
    //def te_opt = te_f != 'NO_FILE' ? te_f : 'NO_FILE'
    //def gene_opt = gene_f != 'NO_FILE' ? gene_f : 'NO_FILE'
    """
    Rscript ${workflow.projectDir}/bin/score_centromeric_classes.R \
        ${repeats_r} ${arrays_r} ${genome_classes} \
        ${metadata} ${te_f} ${gene_f} \
        ${assembly.baseName}_centromeric_scores.csv
    """
}

// ────────────────────── PREDICT CENTROMERIC ──────────────────────
process PREDICT_CENTROMERIC {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path scores
    output:
    path "${scores.baseName}_predictions.csv"
    script:
    """
    python3 ${workflow.projectDir}/bin/predict_centromeric.py \
        ${scores} ${workflow.projectDir}/model/centromeric_model_v2.pkl \
        ${scores.baseName}_predictions.csv
    """
}

// ────────────────────── GC & CTW (parallel) ──────────────────────
process GC {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path assembly
    output:
    path "${assembly.baseName}_GC.csv"
    script:
    """
    Rscript ${workflow.projectDir}/bin/GC.R ${assembly} ${assembly.baseName}_GC.csv
    """
}

process CTW {
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path assembly
    output:
    path "${assembly.baseName}_CTW.csv"
    script:
    """
    Rscript ${workflow.projectDir}/bin/CTW.R ${assembly} ${assembly.baseName}_CTW.csv
    """
}

// ────────────────────── FINAL CAP ──────────────────────
process CAP {
    publishDir "${params.outdir}", mode: 'copy'
	
    input:
    path assembly
    path predictions
    tuple path(repeats_r), path(arrays_r), path(genome_classes)
    path metadata
    val te_f // optional
    val gene_f // optional
    path gc_ch
    path ctw_ch
    path scores
	
    output:
    path "${assembly.baseName}_CAP_plot_*.png"
    path "${assembly.baseName}_CAP_dotplot.png"
    path "${assembly.baseName}_CAP_repeat_families.csv"
    path "${assembly.baseName}_CAP_model.txt"
	
    script:
    """
    Rscript ${workflow.projectDir}/bin/CAP.R \
        ${predictions} ${repeats_r} ${arrays_r} ${genome_classes} \
        ${metadata} ${assembly.baseName} \
        ${gc_ch} ${ctw_ch} ${te_f} ${gene_f} ${scores}
    """
}


// ────────────────────── WORKFLOW ──────────────────────
workflow {
	
// ────────────────────── HELP MESSAGE ──────────────────────
if (!params.assembly) {
    def help = """
╔══════════════════════════════════════════════════════════════════════════════╗
║ CAP - Centromere Analysis Pipeline ║
╚══════════════════════════════════════════════════════════════════════════════╝
Required:
  --assembly <path> Genome assembly in FASTA format
Optional:
  --te_gff <path> EDTA TE annotation (GFF3)
  --gene_gff <path> Helixer gene annotation (GFF)
  --templates <path> Custom repeat templates (FASTA)
  --cores <int> Number of CPU cores [default: ${params.cores}]
  --outdir <path> Output directory [default: ${params.outdir}]
────────────────────────────────────────────────────────────────────────────────
                             INSTALLATION & RUN OPTIONS
────────────────────────────────────────────────────────────────────────────────
1. DOCKER (Recommended – zero setup)
   docker run --rm -v \$(pwd)/data:/data yourname/cap-pipeline:latest \\
     nextflow run main.nf -profile docker --assembly /data/genome.fasta
2. CONDA (Local Linux/macOS/WSL)
   conda env create -f environment.yml
   conda activate cap-pipeline
   nextflow run main.nf -profile conda --assembly data/genome.fasta
3. R + renv (Lightweight, R-only)
   R -e 'renv::restore()' # installs exact R packages
   nextflow run main.nf -profile renv --assembly data/genome.fasta
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

	ready_ch = CHECK_DEPS()
	// Input channels
    assembly_ch = channel.fromPath(params.assembly, checkIfExists: true)
    //templates_ch = params.templates ? channel.fromPath(params.templates) : channel.value(null)
    te_gff_ch = params.te_gff ? channel.fromPath(params.te_gff) : channel.value(null)
    gene_gff_ch = params.gene_gff ? channel.fromPath(params.gene_gff) : channel.value(null)

	// ---- TRASH2 chain ----
	trash_in = assembly_ch.map { assembly ->
    def t = params.templates ? file(params.templates) : null
    tuple(assembly, t)
	}
	trash_out = TRASH2(trash_in, ready_ch)
    filtered = FILTER_TRASH(trash_out)
    reclassed = MERGE_CLASSES(filtered)

    metadata = GET_METADATA(assembly_ch)

    // ---- Optional TEs ----
    if (params.te_gff) {
        te_parsed = PARSE_TES(te_gff_ch)
        te_filtered = FILTER_TES(te_parsed)
    } else {
        te_filtered = channel.of('NO_FILE')
    }

    // ---- Optional Genes ----
    if (params.gene_gff) {
        gene_parsed = PARSE_GENES(gene_gff_ch)
        gene_filtered = FILTER_GENES(gene_parsed)
    } else {
        gene_filtered = channel.of('NO_FILE')
    }

    // ---- Scoring + prediction ----
    scores = SCORE_CENTROMERIC(assembly_ch, reclassed, metadata, te_filtered, gene_filtered)
    preds = PREDICT_CENTROMERIC(scores)

    // ---- GC / CTW in parallel ----
    gc = GC(assembly_ch)
    ctw = CTW(assembly_ch)

    // ---- Final visualisation ----
    CAP(assembly_ch, preds, reclassed, metadata, te_filtered, gene_filtered, gc, ctw, scores)
}