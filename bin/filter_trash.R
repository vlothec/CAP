#!/usr/bin/env Rscript

# filter_trash.R
# Description: Filter TRASH2 outputs

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: filter_trash.R <repeats_with_seq.csv> <arrays.csv> <output_repeats_filtered.csv> <output_arrays_filtered.csv>")
}

repeats_with_seq_csv <- args[1]
arrays_csv <- args[2]
output_repeats_filtered <- args[3]
output_arrays_filtered <- args[4]

# Load data
repeats_data <- read.csv(repeats_with_seq_csv)
arrays_data <- read.csv(arrays_csv)

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Filtering
{
  repeats = repeats_data
  
  arrays = arrays_data
  
  arrays$arrayID = seq_len(nrow(arrays))
  
  
  repeats$arrayID = 0
  for(j in 1 : nrow(arrays)) {
    repeats$arrayID[repeats$seqID == arrays$seqID[j] & 
                      repeats$start >= arrays$start[j] & 
                      repeats$end <= arrays$end[j]] = j
  }
  
  ### Remove repeats (with no class) with score over 60 =====
  nrow(repeats) # 139458
  maxed = 60
  repeats = repeats[repeats$score <= maxed, ]
  
  ### remove arrays and repeats with less than 1000 bp of repeats =========
  min_rep_bp_array <- 1000
  nrow(repeats) # 136956
  nrow(arrays)  # 17937
  
  cat('C:\n')
  arrays_repeat_bp <- unlist(lapply(1 : nrow(arrays), function(X) {
    sum(repeats$width[repeats$arrayID == arrays$arrayID[X]])
  }   ))
  arrayID_to_remove <- which(arrays_repeat_bp < min_rep_bp_array)
  
  arrays <- arrays[-which(arrays$arrayID %in% arrayID_to_remove),]
  repeats <- repeats[-which(repeats$arrayID %in% arrayID_to_remove),]
  
  ### remove arrays and repeats with less than 15 repeats, when average adist score is higher than 25
  nrow(repeats) # 119060  9976
  nrow(arrays)  # 1321   682
  
  cat('D:\n')
  minr = 15
  maxed = 25
  
  arrays_repeat_info <- lapply(seq_len(nrow(arrays)), function(x) {
    which_repeats <- which(repeats$arrayID == arrays$arrayID[x])
    list(
      count = length(which_repeats),
      mean_score = if (length(which_repeats) > 0) mean(repeats$score[which_repeats]) else NA
    )
  })
  repeat_counts <- sapply(arrays_repeat_info, `[[`, "count")
  repeat_means  <- sapply(arrays_repeat_info, `[[`, "mean_score")
  arrayID_to_remove <- arrays$arrayID[repeat_counts < minr & repeat_means > maxed]
  
  arrays  <- arrays[!(arrays$arrayID %in% arrayID_to_remove), ]
  repeats <- repeats[!(repeats$arrayID %in% arrayID_to_remove), ]
  
  nrow(repeats) #  9332
  nrow(arrays)  #  544
  
  
  
}

# Save outputs
write.csv(repeats, file = output_repeats_filtered, row.names = FALSE)
write.csv(arrays, file = output_arrays_filtered, row.names = FALSE)