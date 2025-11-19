#!/usr/bin/env Rscript

# filter_TEs.R
# Description: Filter parsed TEs

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: filter_TEs.R <TEs_parsed.csv> <arrays_filtered.csv> <metadata_file.csv> <output_TEs_filtered.csv>")
}

tes_parsed_csv <- args[1]
arrays_filtered_csv <- args[2]
metadata_file_csv <- args[3]
output_tes_filtered <- args[4]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)
                  library(GenomicRanges)
                  library(IRanges)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# variables
tes_overlap_max_perc = 0.8

# Load data 
tes_data <- read.csv(tes_parsed_csv)
arrays_filtered <- read.csv(arrays_filtered_csv)

metadata <- read.csv(metadata_file_csv)

chromosomes <- metadata$chromosome.name
chromosomes_lengths <- metadata$size

tes_filtered_total = data.frame()
arrays_filtered$overlapping_bp = 0
tes_data$overlapping_bp = 0
# for each chromosome
for (j in seq_along(chromosomes)) {  
  sequence_arrays = arrays_filtered[arrays_filtered$seqID == chromosomes[j], ]
  sequence_tes = tes_data[tes_data$V1 == chromosomes[j], ]
  
  if(nrow(sequence_arrays) + nrow(sequence_tes) == 0) next
  
  if(nrow(sequence_arrays) & nrow(sequence_tes)) {
    gr1 <- with(sequence_arrays, GRanges(chromosomes[j], IRanges(start, end)))
    gr2 <- with(sequence_tes, GRanges(chromosomes[j], IRanges(start, end)))
    
    overlaps <- as.data.frame(findOverlaps(gr1, gr2))
    
    pb <- txtProgressBar(min = 0, max = nrow(overlaps), style = 3)
    for(k in seq_len(nrow(overlaps))) {
      overlap_bp = width(pintersect(gr1[overlaps$queryHits[k]], gr2[overlaps$subjectHits[k]]))
      sequence_arrays$overlapping_bp[overlaps$queryHits[k]] = sequence_arrays$overlapping_bp[overlaps$queryHits[k]] + overlap_bp
      sequence_tes$overlapping_bp[overlaps$subjectHits[k]] = sequence_tes$overlapping_bp[overlaps$subjectHits[k]] + overlap_bp
      setTxtProgressBar(pb, k)
    }
    close(pb)
    sequence_arrays$width = sequence_arrays$end - sequence_arrays$start + 1
    sequence_tes$width = sequence_tes$end - sequence_tes$start + 1
    sequence_tes$overlapping_percentage = sequence_tes$overlapping_bp / sequence_tes$width
    sequence_arrays$overlapping_percentage = sequence_arrays$overlapping_bp / sequence_arrays$width
    
    tes_filtered = sequence_tes[sequence_tes$overlapping_percentage <= tes_overlap_max_perc,]
    
    tes_filtered_total = rbind(tes_filtered_total, tes_filtered)
    
  } else {
    sequence_tes$width = sequence_tes$end - sequence_tes$start + 1
    sequence_tes$overlapping_percentage = rep(0, nrow(sequence_tes))
    tes_filtered = sequence_tes
    tes_filtered_total = rbind(tes_filtered_total, tes_filtered)
    
  }
}

# Save output
write.csv(tes_filtered_total, file = output_tes_filtered, row.names = FALSE)