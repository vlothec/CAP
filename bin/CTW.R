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
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data 
fasta <- read.fasta(assembly_fasta)
assembly.name <- basename(assembly_fasta)

metadata_data <- data.frame(assembly.name = rep(assembly.name, length(fasta)),
                            chromosome.name = names(fasta),
                            size = unlist(lapply(fasta, length)),
                            is.chr = rep(1, length(fasta)))

ctw_data <- data.frame(chromosome = vector(mode = "character"),
                       bin_mid = vector(mode = "numeric"),
                       bin_value = vector(mode = "numeric"))

bin_edta <- 100000 # should be the same as in CAP.R 

for(i in seq_along(fasta)) {
  len <- metadata_data$size[i]
  
  
  bins <- round(len / bin_edta)
  bin_size <- len / bins
  window.starts.CTW = round((0 : bins) * bin_size)
  window.ends.CTW <- window.starts.CTW[-1]
  window.starts.CTW <- c(window.starts.CTW[-length(window.starts.CTW)]) + 1
  CTW.values =  lapply(seq_along(window.starts.CTW), function(X) CTW(paste(fasta[[i]][window.starts.CTW[X] : window.ends.CTW[X]], collapse = ""), depth = 10, desired_alphabet = NULL, beta = NULL))
  log2e <- log(2,base=exp(1))
  CTW.values <- unlist(CTW.values)
  CTW.values = CTW.values * -log2e
  CTW.values <- CTW.values / (bin_edta - 10)
  CTW.values[CTW.values > 1] = 1
  CTW.values[CTW.values < 0] = 0
  CTW.values = CTW.values * 100
  
  mids <- window.starts.CTW + bin_edta/2
  
  ctw_data <- rbind(ctw_data, data.frame(chromosome = rep(metadata_data$chromosome.name[i], length(mids)),
                                         bin_mid = mids,
                                         bin_value = CTW.values))
  
}


# --- Your CTW computation code here ---

# Save output
write.csv(ctw_data, file = output_ctw, row.names = FALSE)  # Assuming ctw_data is created in your code