#!/usr/bin/env Rscript

# CTW.R
# Description: Compute CTW from assembly (assuming some metric)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: CTW.R <assembly.fasta> <output_CTW.csv>")
}

assembly_fasta <- args[1]
output_ctw <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source("./aux_fun.R")

# Load data 


# --- Your CTW computation code here ---

# Save output
write.csv(ctw_data, file = output_ctw, row.names = FALSE)  # Assuming ctw_data is created in your code