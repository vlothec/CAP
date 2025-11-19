#!/usr/bin/env Rscript

# filter_genes.R
# Description: Filter parsed genes

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: filter_genes.R <genes_parsed.csv> <output_genes_filtered.csv>")
}

genes_parsed_csv <- args[1]
output_genes_filtered <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))
genes_data <- data.frame()

# Load data 

# --- Your filtering code here ---

# Save output
write.csv(genes_data, file = output_genes_filtered, row.names = FALSE)