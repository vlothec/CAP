#!/usr/bin/env Rscript

# filter_genes.R
# Description: Filter parsed genes

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: filter_genes.R <genes_parsed.csv> <arrays_filtered.csv> <metadata_file.csv> <output_genes_filtered.csv>")
}

genes_parsed_csv <- args[1]
arrays_filtered_csv <- args[2]
metadata_file_csv <- args[3]
output_genes_filtered <- args[4]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)
                  library(GenomicRanges)
                  library(IRanges)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# variables
genes_overlap_max_perc <- 0.8

# Load data
genes_data <- read.csv(genes_parsed_csv)
arrays_filtered <- read.csv(arrays_filtered_csv)
metadata <- read.csv(metadata_file_csv)

chromosomes <- metadata$chromosome.name
chromosomes_lengths <- metadata$size

genes_filtered_total <- data.frame()
arrays_filtered$overlapping_bp <- 0
genes_data$overlapping_bp <- 0
genes_data$width <- genes_data$end - genes_data$start + 1

for (j in seq_along(chromosomes)) {  
  sequence_arrays = arrays_filtered[arrays_filtered$seqID == chromosomes[j], ]
  sequence_genes <- genes_data[genes_data$seqID == chromosomes[j],]
  
  if(nrow(sequence_genes) == 0) next
  if(nrow(sequence_arrays) == 0) {
    genes_filtered_total = rbind(genes_filtered_total, sequence_genes)
    next
  }
  
  gene_IDs <- unique(sequence_genes$ID[sequence_genes$V3 == "gene"])
  gene_cluster_start <- which(sequence_genes$V3 == "gene")
  gene_cluster_end = c(gene_cluster_start[-1], nrow(sequence_genes))
  
  which_sequence_CDS <- which(sequence_genes$V3 == "CDS")
  sequence_genes_CDS = sequence_genes[which_sequence_CDS, ]
  sequence_genes_CDS$overlapping_bp <- 0
  
  gr1 <- with(sequence_arrays, GRanges( chromosomes[j], IRanges(start, end)))
  gr2 <- with(sequence_genes_CDS, GRanges(chromosomes[j], IRanges(V4, V5)))
  
  overlaps <- as.data.frame(findOverlaps(gr1, gr2))
  
  if (nrow(overlaps) > 0) {
    overlap_bp_vec <- width(pintersect(gr1[overlaps$queryHits], gr2[overlaps$subjectHits]))
    overlap_sums <- aggregate(overlap_bp_vec, by = list(subjectHits = overlaps$subjectHits), FUN = sum)
    sequence_genes_CDS$overlapping_bp[overlap_sums$subjectHits] <- sequence_genes_CDS$overlapping_bp[overlap_sums$subjectHits] + overlap_sums$x
  }
  sequence_genes$overlapping_bp[which_sequence_CDS] <- sequence_genes_CDS$overlapping_bp
  genes_to_remove <- NULL
  
  for(k in seq_along(gene_IDs)) { 
    gene_cluster <- sequence_genes[gene_cluster_start[k] : gene_cluster_end[k], ]
    gene_CDS_total_size <- sum(gene_cluster$width[gene_cluster$V3 == "CDS"])
    gene_CDS_total_overlap <- sum(gene_cluster$overlapping_bp)
    if((gene_CDS_total_overlap / gene_CDS_total_size) > genes_overlap_max_perc) {
      genes_to_remove <- c(genes_to_remove, (gene_cluster_start[k] : gene_cluster_end[k]))
    }
  }
  if(!is.null(genes_to_remove)) {
    genes_filtered <- sequence_genes[-genes_to_remove, ]
  } else {
    genes_filtered <- sequence_genes
  }
  
  
  genes_filtered_total = rbind(genes_filtered_total, genes_filtered)
  
  
}

# Save output
write.csv(genes_filtered_total, file = output_genes_filtered, row.names = FALSE)