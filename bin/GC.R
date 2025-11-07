#!/usr/bin/env Rscript

# GC.R
# Description: Compute GC content from assembly

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: GC.R <assembly.fasta> <output_GC.csv>")
}

assembly_fasta <- args[1]
output_gc <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source("./aux_fun.R")

# Load data

# --- Your GC computation code here ---

# Save output
write.csv(gc_data, file = output_gc, row.names = FALSE)  # Assuming gc_data is created in your code