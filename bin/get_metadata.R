#!/usr/bin/env Rscript

# get_metadata.R
# Description: Extract metadata from assembly

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: get_metadata.R <assembly.fasta> <output_metadata.csv>")
}

assembly_fasta <- args[1]
output_metadata <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data 

fasta <- read.fasta(assembly_fasta)

# Run
{
  assembly.name <- basename(assembly_fasta)
  
  metadata_data <- data.frame(assembly.name = rep(assembly.name, length(fasta)),
                              chromosome.name = names(fasta),
                              size = unlist(lapply(fasta, length)),
                              is.chr = rep(1, length(fasta)))
}

# Save output
write.csv(metadata_data, file = output_metadata, row.names = FALSE) 