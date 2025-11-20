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
suppressMessages({library(Biostrings)
                  library(BCT)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data 
fasta <- readDNAStringSet(assembly_fasta)
assembly.name <- basename(assembly_fasta)

metadata_data <- data.frame(assembly.name = rep(assembly.name, length(fasta)),
                            chromosome.name = names(fasta),
                            size = width(fasta),
                            is.chr = rep(1, length(fasta)))
for(i in seq_len(nrow(metadata_data))) {
  metadata_data$chromosome.name[i] <- strsplit(metadata_data$chromosome.name[i], " ")[[1]][1]
}
ctw_data <- data.frame(chromosome = vector(mode = "character"),
                       bin_mid = vector(mode = "numeric"),
                       bin_value = vector(mode = "numeric"))

bin_edta <- 100000 # should be the same as in CAP.R 

for(i in seq_along(fasta)) {
  print(paste0("Processing chromosome: ", metadata_data$chromosome.name[i]))
  len <- metadata_data$size[i]
  
  
  bins <- max(1, round(len / bin_edta))
  bin_size <- len / bins
  window.starts.CTW = round((0 : bins) * bin_size)
  window.ends.CTW <- window.starts.CTW[-1]
  window.starts.CTW <- c(window.starts.CTW[-length(window.starts.CTW)]) + 1
  CTW.values =  lapply(seq_along(window.starts.CTW), function(X) CTW(as.character(subseq(fasta[[i]], window.starts.CTW[X], window.ends.CTW[X])), depth = 10, desired_alphabet = NULL, beta = NULL))
  log2e <- log(2,base=exp(1))
  CTW.values <- unlist(CTW.values)
  CTW.values = CTW.values * -log2e
  actual_bin_sizes <- window.ends.CTW - window.starts.CTW + 1
  CTW.values <- CTW.values / (actual_bin_sizes - 10)
  CTW.values[CTW.values > 1] = 1
  CTW.values[CTW.values < 0] = 0
  CTW.values = CTW.values * 100
  mids <- (window.starts.CTW + window.ends.CTW) / 2
  mids <- window.starts.CTW + bin_edta/2
  
  ctw_data <- rbind(ctw_data, data.frame(chromosome = rep(metadata_data$chromosome.name[i], length(mids)),
                                         bin_mid = mids,
                                         bin_value = CTW.values))
  
}


# Save output
write.csv(ctw_data, file = output_ctw, row.names = FALSE)