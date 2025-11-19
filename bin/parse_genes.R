#!/usr/bin/env Rscript

# parse_genes.R
# Description: Parse gene annotation from GFF

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: parse_genes.R <helixer.gff> <output_genes_parsed.csv>")
}

gene_gff <- args[1]
output_genes_parsed <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data 

# --- Your parsing code here ---
parsed_data <- data.frame()

# Save output
write.csv(parsed_data, file = output_genes_parsed, row.names = FALSE)  # Assuming parsed_data is created in your code