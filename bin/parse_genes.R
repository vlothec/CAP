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
gene_raw_data <- read.table(gene_gff,
                            header = FALSE,
                            sep = "\t",
                            comment.char = "#",
                            blank.lines.skip = TRUE,
                            stringsAsFactors = FALSE)

# --- Your parsing code here ---
names(gene_raw_data) <- c("seqID", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

attributes = lapply(gene_raw_data$attributes, function(X) strsplit(X, split = ";")[[1]])
new_cols = lapply(attributes, function(X) unlist(lapply(X, function(x) strsplit(x, split = "=")[[1]][1])))
new_cols = unique(unlist(new_cols))
gene_full = gene_raw_data[, 1:9]
for (j in seq_along(new_cols)) {
  cat(j, "/", length(new_cols), "\n")
  new_data = lapply(gene_raw_data$attributes, function(X) {
  new_data = unlist(lapply(new_data, function(X) {
    split_res <- strsplit(X, split = ";")[[1]]
    if (length(split_res) > 0) split_res[1] else NA
  }))
    if (length(m) > 0) sub(paste0(new_cols[j], "="), "", m) else NA
  })
  new_data = unlist(new_data)
  
  gene_full = cbind(gene_full, new_data)
  names(gene_full)[ncol(gene_full)] = new_cols[j]
}



# Save output
write.csv(parsed_data, file = output_genes_parsed, row.names = FALSE)  # Assuming parsed_data is created in your code