#!/usr/bin/env Rscript

# merge_classes.R
# Description: Merge and reclassify filtered repeats and arrays

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: merge_classes.R <repeats_filtered.csv> <arrays_filtered.csv> <output_repeats_reclassed.csv> <output_arrays_reclassed.csv> <output_genome_classes.csv>")
}

repeats_filtered_csv <- args[1]
arrays_filtered_csv <- args[2]
output_repeats_reclassed <- args[3]
output_arrays_reclassed <- args[4]
output_genome_classes <- args[5]

# Load data
repeats_data <- read.csv(repeats_filtered_csv)
arrays_data <- read.csv(arrays_filtered_csv)

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Merging
{
  #TODO update
  repeats <- repeats_data
  arrays <- arrays_data
  

  classes = data.frame(class = unique(repeats$class))
  
  classes$count = 0
  classes$consensus = ""
  classes$median_length = 0
  classes$length_SD = 0
  classes$mean_edit_score = 0 
  classes$total_bp = 0
  
  min_repeats_to_align = 10
  
  for(j in seq_len(nrow(classes))) {
    repeats_class = repeats[repeats$class == classes$class[j], ]
    
    class_name_nchar = as.numeric(strsplit(classes$class[j], split = "_")[[1]][1])
    
    classes$count[j] = nrow(repeats_class)
    classes$median_length[j] = round(median(repeats_class$width))
    classes$length_SD[j] = sd(repeats_class$width)
    classes$total_bp[j] = sum(repeats_class$width)
    
    if(nrow(repeats_class) < min_repeats_to_align) {
      repeats_to_align = seq_len(nrow(repeats_class))
    } else if(nrow(repeats_class) < 10000) {
      repeats_to_align = sample(seq_len(nrow(repeats_class)), (nrow(repeats_class) %/% 100 + min_repeats_to_align), replace = FALSE)
    } else if(nrow(repeats_class) < 100000) {
      repeats_to_align = sample(seq_len(nrow(repeats_class)), (nrow(repeats_class) %/% 1000 + min_repeats_to_align), replace = FALSE)
    } else {
      repeats_to_align = sample(seq_len(nrow(repeats_class)), (nrow(repeats_class) %/% 10000 + min_repeats_to_align), replace = FALSE)
    }
    
    sequences_to_align <- repeats_class$sequence[repeats_to_align]
    if(length(sequences_to_align) == 0) {
      cat("\n\n\n\n\n\n\n\n\n\n")
      stop(paste0(" did not find repeats in one of the classes: ", classes$class[j], ", investigate"))
    } else if(length(sequences_to_align) == 1) {
      classes$consensus[j] <- sequences_to_align
      classes$mean_edit_score[j] = 0
    } else {
      a <- capture.output({alignment_matrix = msa(sequences_to_align, method = "ClustalOmega", type = "dna")})
      classes$consensus[j] <- consensus_N(alignment_matrix, class_name_nchar)
      classes$mean_edit_score[j] = mean(adist(classes$consensus[j], sequences_to_align))
    }
  }
  
  ### save classes
  classes$importance = classes$count * classes$median_length
  
  classes$new_class = ""
  classes$new_importance = 0
  
  ### remove classes with less than 5000 bp of repeats for merging, to speed it up ========
  if(nrow(classes[classes$total_bp >= 5000, ]) >= 20) {
    classes_discarded = classes[classes$total_bp < 5000, ]
    classes = classes[classes$total_bp >= 5000, ]
    if(nrow(classes_discarded) > 0) {
      classes_discarded$new_class = classes_discarded$class
      classes_discarded$score = 0
      classes_discarded$sum_coverage = classes_discarded$importance
    }
  } else if(nrow(classes[classes$total_bp >= 1000, ]) >= 20) {
    classes_discarded = classes[classes$total_bp < 1000, ]
    classes = classes[classes$total_bp >= 1000, ]
    if(nrow(classes_discarded) > 0) {
      classes_discarded$new_class = classes_discarded$class
      classes_discarded$score = 0
      classes_discarded$sum_coverage = classes_discarded$importance
    }
  } else if(nrow(classes[classes$total_bp >= 500, ]) >= 20) {
    classes_discarded = classes[classes$total_bp < 1000, ]
    classes = classes[classes$total_bp >= 500, ]
    if(nrow(classes_discarded) > 0) {
      classes_discarded$new_class = classes_discarded$class
      classes_discarded$score = 0
      classes_discarded$sum_coverage = classes_discarded$importance
    }
  } else {
    classes_discarded = classes[NULL, ]
  }
  
  ### check kmer based similarity between all ==============
  
  
  classes = classes[order(classes$importance, decreasing = TRUE), ]
  
  similarity_matrix = matrix(nrow = nrow(classes), ncol = nrow(classes))
  for(j in 1 : nrow(similarity_matrix)) {
    for(k in j : ncol(similarity_matrix)) {
      similarity_matrix[j,k] = kmer_compare(classes$consensus[j], classes$consensus[k])
      similarity_matrix[k,j] = kmer_compare(classes$consensus[k], classes$consensus[j])
    }
  }
  
  
  for(j in 1 : nrow(classes)) {
    if(classes$new_class[j] == "") {
      classes$new_class[j] = classes$class[j]
      other_similar = which(similarity_matrix[j, ] > 0.2 & similarity_matrix[, j] > 0.2)
      classes$new_class[other_similar] = classes$class[j]
      classes$new_importance[j] = sum(classes$importance[other_similar])
    }
  }
  
  classes_new = classes[classes$new_importance != 0, ]
  classes_new$importance = classes_new$new_importance
  classes_new$class = classes_new$new_class
  classes_new$score = classes_new$mean_edit_score / classes_new$median_length
  classes_new$sum_coverage = classes_new$importance
  
  if(nrow(classes_discarded) != 0) {
    classes_new = rbind(classes_new, classes_discarded)
  } 
  
  ### Assign new classes to repeats and arrays =============
  
  classes_new$class_num_ID = 1 : nrow(classes_new)
  classes$class_num_ID = 0
  for(j in 1 : nrow(classes)) {
    classes$class_num_ID[j] = classes_new$class_num_ID[which(classes_new$new_class == classes$new_class[j])]
  }
  
  repeats$new_class = ""
  repeats$new_class_num_ID = 0
  arrays$new_class_num_ID = 0
  for(j in 1 : nrow(classes)){
    repeats$new_class[repeats$class == classes$class[j]] = classes$new_class[j]
    repeats$new_class_num_ID[repeats$class == classes$class[j]] = classes$class_num_ID[j]
    arrays$new_class_num_ID[arrays$class == classes$class[j]] = classes$class_num_ID[j]
  }
  repeats$new_class[repeats$new_class == ""] = repeats$class[repeats$new_class == ""]
  
  ### Shrink arrays to only reach over repeats =============
  
  
  for(j in seq_len(nrow(arrays))) {
    repeats_array = repeats[repeats$arrayID == arrays$arrayID[j], ]
    arrays$start[j] = min(repeats_array$start)
    arrays$end[j] = max(repeats_array$end)
  }
  
  
  
  
}

# Save outputs
write.csv(repeats, file = output_repeats_reclassed, row.names = FALSE)
write.csv(arrays, file = output_arrays_reclassed, row.names = FALSE)
write.csv(classes_new, file = output_genome_classes, row.names = FALSE)  