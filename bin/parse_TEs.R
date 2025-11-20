#!/usr/bin/env Rscript

# parse_TEs.R
# Description: Parse TE annotation from GFF3 and output as CSV

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: parse_TEs.R <te_annotation.gff3> <metadata.csv> <output_TEs_parsed.csv> \nNote: Output will be written in CSV format.")
}

te_gff <- args[1]
metadata_csv <- args[2]
output_tes_parsed <- args[3]
print(te_gff)
print(metadata_csv)
print(output_tes_parsed)
if (!file.exists(metadata_csv)) {
  stop(paste("Metadata CSV file not found:", metadata_csv))
}
if (!file.exists(te_gff)) {
  stop(paste("TE GFF file not found:", te_gff))
}

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
metadata <- read.csv(metadata_csv)

# Rename TE seqID if needed
metadata$new_chr_name <- ""
metadata$shortened_chr_name = ""
for (i in seq_along(unique(metadata$chromosome))) {
  metadata$new_chr_name[metadata$chromosome == unique(metadata$chromosome)[i]] <- paste0("Chr", i)
  metadata$shortened_chr_name[metadata$chromosome == unique(metadata$chromosome)[i]] <- strsplit(unique(metadata$chromosome)[i], split = " ")[[1]][1]
  te_raw_data$seqID[te_raw_data$seqID == paste0("Chr", i)] <- strsplit(unique(metadata$chromosome)[i], split = " ")[[1]][1]
}

# Parse attributes and reformat data
attributes = lapply(te_raw_data$attributes, function(X) strsplit(X, split = ";")[[1]])
new_cols = lapply(attributes, function(X) unlist(lapply(X, function(x) strsplit(x, split = "=")[[1]][1])))
new_cols = unique(unlist(new_cols))
edta_full = te_raw_data[, 1:9]
for (j in seq_along(new_cols)) {
  new_data <- sapply(te_raw_data$attributes, function(X) {
    split_res <- strsplit(X, split = ";")[[1]]
    m <- split_res[grep(paste0("^", new_cols[j], "="), split_res)]
    if (length(m) > 0) sub(paste0(new_cols[j], "="), "", m[1]) else NA
  })
  
  edta_full <- cbind(edta_full, new_data)
  names(edta_full)[ncol(edta_full)] <- new_cols[j]
}
edta_full$old_type = edta_full$type

names(edta_full) <- tolower(names(edta_full))
names(edta_full)[names(edta_full) == "seqid"] <- "seqID"

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