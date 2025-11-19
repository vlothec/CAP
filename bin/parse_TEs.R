#!/usr/bin/env Rscript

# parse_TEs.R
# Description: Parse TE annotation from GFF3

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: parse_TEs.R <EDTA.TEanno.split.gff3> <output_TEs_parsed.csv>")
}

te_gff <- args[1]
output_tes_parsed <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data 
parsed_data <- data.frame()


# --- Your parsing code here ---

# Save output
write.csv(parsed_data, file = output_tes_parsed, row.names = FALSE)  # Assuming parsed_data is created in your code