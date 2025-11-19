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
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data
fasta <- read.fasta(assembly_fasta)
assembly.name <- basename(assembly_fasta)

metadata_data <- data.frame(assembly.name = rep(assembly.name, length(fasta)),
                            chromosome.name = names(fasta),
                            size = unlist(lapply(fasta, length)),
                            is.chr = rep(1, length(fasta)))

gc_data <- data.frame(chromosome = vector(mode = "character"),
                       bin_mid = vector(mode = "numeric"),
                       bin_value = vector(mode = "numeric"))

bin_gc <- 2000 # should be the same as in CAP.R 

for(i in seq_along(fasta)) {
  len <- metadata_data$size[i]
  win_gc <- genomic.bins.starts(1, len, bin.size = bin_gc)
  gc_vals <- calculate.GC.in.windows.2(win_gc, fasta[[i]], bin_gc)
  mids <- win_gc + bin_gc/2
  
  gc_data <- rbind(gc_data, data.frame(chromosome = rep(metadata_data$chromosome.name[i], length(mids)),
                                       bin_mid = mids,
                                       bin_value = gc_vals))
  
}

# Save output
write.csv(gc_data, file = output_gc, row.names = FALSE)