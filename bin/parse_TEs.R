#!/usr/bin/env Rscript

# parse_TEs.R
# Description: Parse TE annotation from GFF3 and output as CSV

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: parse_TEs.R <te_annotation.gff3> <output_TEs_parsed.csv>\nNote: Output will be written in CSV format.")
}

te_gff <- args[1]
output_tes_parsed <- args[2]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load data
te_raw_data <- read.table(te_gff,
                          header = FALSE,
                          sep = "\t",
                          comment.char = "#", 
                          blank.lines.skip = TRUE,
                          stringsAsFactors = FALSE)
names(te_raw_data) <- c("seqID", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

attributes = lapply(te_raw_data$attributes, function(X) strsplit(X, split = ";")[[1]])
new_cols = lapply(attributes, function(X) unlist(lapply(X, function(x) strsplit(x, split = "=")[[1]][1])))
new_cols = unique(unlist(new_cols))
edta_full = te_raw_data[, 1:9]
for (j in seq_along(new_cols)) {
  cat(j, "/", length(new_cols), "\n")
  new_data = lapply(te_raw_data$attributes, function(X) {
  new_data = unlist(lapply(new_data, function(X) {
    split_res <- strsplit(X, split = ";")[[1]]
    if (length(split_res) > 0) split_res[1] else NA
  }))
    if (length(m) > 0) sub(paste0(new_cols[j], "="), "", m) else NA
  })
  new_data = unlist(new_data)
  
  edta_full = cbind(edta_full, new_data)
  names(edta_full)[ncol(edta_full)] = new_cols[j]
}
edta_full$old_type = edta_full$type

names(edta_full) <- tolower(names(edta_full))

edta_repeat_region = edta_full[edta_full$type == "repeat_region", ]
edta_non_rep_reg = edta_full[edta_full$type != "repeat_region", ]

names = unique(edta_repeat_region$name)
names = sort(names)

names = data.frame(names)
names$short = unlist(lapply(names$names, function(X) {
  pos = grep("[^A-Za-z]", strsplit(X, split = "")[[1]])
  if(length(pos) == 0 ) return(X)
  return(substr(X, 1, min(pos) - 1))
  } ))

names_short = unique(names$short)

edta_repeat_region$old_type = edta_repeat_region$type
idx <- match(edta_repeat_region$name, names$names)
edta_repeat_region$type <- ifelse(is.na(idx), edta_repeat_region$type, names$short[idx])

edta_modified = rbind(edta_non_rep_reg, edta_repeat_region)

edta_modified$score = "."
edta_modified$phase = "."




# Save output
write.csv(edta_modified, file = output_tes_parsed, row.names = FALSE)  