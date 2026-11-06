#!/usr/bin/env Rscript

# merge_classes.R
# Description: Merge and reclassify filtered repeats and arrays
# Optimized for performance

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
repeats <- read.csv(repeats_filtered_csv)
arrays <- read.csv(arrays_filtered_csv)

# Stabilize random repeat subsampling across runs.
set.seed(0)

# Load libraries
suppressMessages({
  library(seqinr)
  library(msa)
})

# Load additional functions
# Robustly find auxfuns.R in the same directory as this script
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)])
script_dir <- dirname(script_name)
if(length(script_dir) == 0 || script_dir == ".") script_dir <- "bin" # Fallback
source(file.path(script_dir, "auxfuns.R"))

# --- 1. Class Statistics & Consensus ---

classes <- data.frame(class = unique(repeats$class))
classes$count <- 0
classes$consensus <- ""
classes$median_length <- 0
classes$length_SD <- 0
classes$mean_edit_score <- 0 
classes$total_bp <- 0

min_repeats_to_align <- 10
max_repeats_to_align <- 100 # Cap to prevent MSA from hanging

for(j in seq_len(nrow(classes))) {
  # Subset repeats for this class
  repeats_class <- repeats[repeats$class == classes$class[j], ]
  
  class_name_nchar <- as.numeric(strsplit(classes$class[j], split = "_")[[1]][1])
  
  classes$count[j] <- nrow(repeats_class)
  classes$median_length[j] <- round(median(repeats_class$width))
  classes$length_SD[j] <- sd(repeats_class$width)
  classes$total_bp[j] <- sum(repeats_class$width)
  
  # Sampling logic (Simplified and Capped)
  n_reps <- nrow(repeats_class)
  if(n_reps <= min_repeats_to_align) {
    repeats_to_align <- seq_len(n_reps)
  } else {
    # Scale sample size but cap at max_repeats_to_align
    sample_size <- min(max_repeats_to_align, (n_reps %/% 100 + min_repeats_to_align))
    repeats_to_align <- sample(seq_len(n_reps), sample_size, replace = FALSE)
  }
  
  sequences_to_align <- repeats_class$sequence[repeats_to_align]
  
  if(length(sequences_to_align) == 0) {
    stop(paste0("Error: No sequences found for class: ", classes$class[j]))
  } else if(length(sequences_to_align) == 1) {
    classes$consensus[j] <- sequences_to_align
    classes$mean_edit_score[j] <- 0
  } else {
    # Capture output to silence MSA messages
    capture.output({
      alignment_matrix <- msa::msa(sequences_to_align, method = "ClustalOmega", type = "dna")
    })
    classes$consensus[j] <- consensus_N(alignment_matrix, class_name_nchar)
    classes$mean_edit_score[j] <- mean(adist(classes$consensus[j], sequences_to_align))
  }
}

classes$importance <- classes$count * classes$median_length
classes$new_class <- ""
classes$new_importance <- 0

# --- 2. Filtering Classes for Merging ---

# Determine threshold based on data distribution
threshold <- 5000
if(nrow(classes[classes$total_bp >= 5000, ]) < 20) {
  if(nrow(classes[classes$total_bp >= 1000, ]) >= 20) {
    threshold <- 1000
  } else {
    threshold <- 500
  }
}

classes_discarded <- classes[classes$total_bp < threshold, ]
classes <- classes[classes$total_bp >= threshold, ]

if(nrow(classes_discarded) > 0) {
  classes_discarded$new_class <- classes_discarded$class
  classes_discarded$score <- 0
  classes_discarded$sum_coverage <- classes_discarded$importance
}

# --- 3. Similarity & Merging ---

classes <- classes[order(classes$importance, decreasing = TRUE), ]

# Calculate similarity matrix
n_classes <- nrow(classes)
similarity_matrix <- matrix(0, nrow = n_classes, ncol = n_classes)

if (n_classes > 0) {
    for(j in 1:n_classes) {
      for(k in j:n_classes) {
        sim <- kmer_compare(classes$consensus[j], classes$consensus[k])
        similarity_matrix[j,k] <- sim
        similarity_matrix[k,j] <- sim
      }
    }
    
    # Greedy merging
    for(j in 1:n_classes) {
      if(classes$new_class[j] == "") {
        # This class becomes a new representative
        classes$new_class[j] <- classes$class[j]
        
        # Find similar classes
        other_similar <- which(similarity_matrix[j, ] > 0.2)
        
        classes$new_class[other_similar] <- classes$class[j]
        classes$new_importance[j] <- sum(classes$importance[other_similar])
      }
    }
}

classes_new <- classes[classes$new_importance != 0, ]
classes_new$importance <- classes_new$new_importance
classes_new$class <- classes_new$new_class
classes_new$score <- classes_new$mean_edit_score / classes_new$median_length
classes_new$sum_coverage <- classes_new$importance

if(nrow(classes_discarded) > 0) {
  classes_new <- rbind(classes_new, classes_discarded)
} 

# --- 4. Map New Classes (Vectorized) ---

# Create mapping vectors
all_mappings <- rbind(
    classes[, c("class", "new_class")],
    classes_discarded[, c("class", "new_class")]
)
map_vec <- setNames(all_mappings$new_class, all_mappings$class)

# Assign IDs
classes_new$class_num_ID <- 1:nrow(classes_new)
id_map_vec <- setNames(classes_new$class_num_ID, classes_new$class)

# Apply mappings
repeats$new_class <- map_vec[repeats$class]
repeats$new_class_num_ID <- id_map_vec[repeats$new_class]
arrays$new_class_num_ID <- id_map_vec[map_vec[arrays$class]]

# Handle NAs if any
repeats$new_class[is.na(repeats$new_class)] <- repeats$class[is.na(repeats$new_class)]

# --- 5. Shrink Arrays (Vectorized) ---

# Calculate min start and max end for each arrayID
starts <- aggregate(start ~ arrayID, data = repeats, min)
ends <- aggregate(end ~ arrayID, data = repeats, max)

# Map back to arrays dataframe
idx <- match(arrays$arrayID, starts$arrayID)
arrays$start <- starts$start[idx]
arrays$end <- ends$end[idx]

# Save outputs
write.csv(repeats, file = output_repeats_reclassed, row.names = FALSE)
write.csv(arrays, file = output_arrays_reclassed, row.names = FALSE)
write.csv(classes_new, file = output_genome_classes, row.names = FALSE)  