#!/usr/bin/env Rscript

# filter_TEs.R
# Description: Filter parsed TEs

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: filter_TEs.R <TEs_parsed.csv> <output_TEs_filtered.csv>")
}

tes_parsed_csv <- args[1]
output_tes_filtered <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))
tes_data <- data.frame()

# Load data 

# --- Your filtering code here ---

# Save output
write.csv(tes_data, file = output_tes_filtered, row.names = FALSE)