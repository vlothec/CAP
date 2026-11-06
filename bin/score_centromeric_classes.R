#!/usr/bin/env Rscript

# score_centromeric_classes.R
# Description: Score centromeric classes with optional TEs and genes

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) > 7) {
  stop("Usage: score_centromeric_classes.R <repeats_reclassed.csv> <arrays_reclassed.csv> <genome_classes.csv> <metadata.csv> [<TEs_filtered.csv>] [<genes_filtered.csv>] <output_centromeric_scores.csv>")
}
no_edta <- FALSE
no_heli <- FALSE

repeats_reclassed_csv <- args[1]
arrays_reclassed_csv <- args[2]
genome_classes_csv <- args[3]
metadata_csv <- args[4]
tes_filtered_csv <- if (args[5] != "NO_FILE") args[5] else no_edta <- TRUE
genes_filtered_csv <- if (args[6] != "NO_FILE") args[6] else no_heli <- TRUE
output_centromeric_scores <- args[7]

# Stabilize random repeat subsampling across runs.
set.seed(0)

# Load libraries
suppressMessages({library(seqinr)
  library(msa)
  library(GenomicRanges)})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load required data
# print("load repeats")
repeats <- read.csv(repeats_reclassed_csv)
# print("load arrays")
arrays <- read.csv(arrays_reclassed_csv)
# print("load classes")
genome_classes_data <- read.csv(genome_classes_csv)
# print("load metadata")
metadata_data <- read.csv(metadata_csv)
for(i in seq_len(nrow(metadata_data))) {
  metadata_data$chromosome.name[i] <- strsplit(metadata_data$chromosome.name[i], " ")[[1]][1]
}

if(nrow(genome_classes_data) == 0) {
  empty_output <- data.frame(
    class = character(), new_class_num_ID = numeric(), count = numeric(),
    mean_length = numeric(), total_bp = numeric(), chromosome = character(),
    total_bp_norm_chr = numeric(), total_bp_norm_rep = numeric(), 
    start_sd_norm_chr = numeric(), start_norm_chr_0_50 = numeric(),
    gaps_count = numeric(), gaps_with_TEs_fraction = numeric(),
    centre_array_edit = numeric(), centre_array_width_sd = numeric(),
    centre_chromosome_edit = numeric(), centre_chromosome_width_sd = numeric(),
    array_sizes_sd_norm_mean_arr_size = numeric(), array_count = numeric(),
    TE_prox_dist = numeric(), TE_prox_SD = numeric(), TE_lm_coef = numeric(),
    TE_prox_score = numeric(), gene_prox_dist = numeric(), gene_prox_SD = numeric(),
    gene_lm_coef = numeric(), gene_prox_score = numeric(),
    pred_centrophilic_TE = character(), t_test_t_val = numeric(), t_test_p_val = numeric()
  )
  write.csv(empty_output, output_centromeric_scores, row.names = FALSE)
  quit(status = 0)
}


# Load optionals if provided
# print("load edta")
if (!no_edta) {
  edta <- read.csv(tes_filtered_csv)
  if (length(edta) == 0 || nrow(edta) == 0) no_edta <- TRUE
}
# print("load genes")
if (!no_heli) {
  genes <- read.csv(genes_filtered_csv)
  if (length(genes) == 0 || nrow(genes) == 0) no_heli <- TRUE
}


# Score
all_scores <- data.frame()
{
  chromosomes <- metadata_data$chromosome.name
  chromosomes_lengths <- metadata_data$size
  
  
  ### Decide which classes to score
  
  # Filtering variables
  MIN_BP_TO_SCORE <- 2000
  MIN_REPEAT_COUNT_TO_SCORE <- 10
  MIN_REP_WIDTH_TO_SCORE <- 8
  N_MAX_REP_PER_CHR <- 4
  MAX_REPEATS_PER_CHROMOSOME <- 8
  
  # find top N_MAX_REP_PER_CHR classes for each chromosome
  top_N_unique <- NULL
  for (j in seq_along(chromosomes)) {  
    sequence_repeats = repeats[repeats$seqID == chromosomes[j], ]
    chromosome_no_rep_size <- chromosomes_lengths[j] - sum(sequence_repeats$width)
    if(nrow(sequence_repeats) == 0) next
    
    chr_classes <- as.data.frame(table(sequence_repeats$new_class))
    names(chr_classes) <- c("class", "count")
    chr_classes$class <- as.character(chr_classes$class)
    chr_classes$mean_length <- 0
    chr_classes$total_bp <- 0
    for(k in seq_len(nrow(chr_classes))) {
      chr_classes$mean_length[k] <- mean(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
      chr_classes$total_bp[k] <- sum(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
    }
    chr_classes <- chr_classes[chr_classes$mean_length >= MIN_REP_WIDTH_TO_SCORE, ]
    chr_classes <- chr_classes[chr_classes$total_bp >= MIN_BP_TO_SCORE, ]
    chr_classes <- chr_classes[chr_classes$count >= MIN_REPEAT_COUNT_TO_SCORE, ]
    if(nrow(chr_classes) == 0) next
    chr_classes <- chr_classes[order(chr_classes$total_bp, decreasing = TRUE), ]
    if(nrow(chr_classes) > 5) top_N_unique <- c(top_N_unique, chr_classes$class[1:5]) else top_N_unique <- c(top_N_unique, chr_classes$class)
    
  }
  top_N_unique <- unique(top_N_unique)
  
  ### Calculate genome_classes characteristics for each repeat class per chromosome
  
  for(j in seq_along(chromosomes)) {
    
    sequence_repeats = repeats[repeats$seqID == chromosomes[j], ]
    sequence_arrays = arrays[arrays$seqID == chromosomes[j], ]
    
    repeats_sorted <- sequence_repeats[order(sequence_repeats$start), ]
    sorted_starts <- repeats_sorted$start
    cum_widths <- cumsum(repeats_sorted$width)
    
    compute_adjustment <- function(pos) {
      idxs <- findInterval(pos, sorted_starts)
      c(0, cum_widths)[idxs + 1]
    }
    
    if(!no_edta) {
      # mask repeats from TE coordinates and find adjusted repeat coordinates
      sequence_edta = edta[edta$seqID == chromosomes[j], ]
      if(nrow(sequence_edta) == 0) next
      gr2 <- with(sequence_edta, GRanges(chromosomes[j], IRanges(start, end)))
      sequence_edta_no_rep <- sequence_edta
      sequence_edta_no_rep$legacy_start <- sequence_edta_no_rep$start
      sequence_edta_no_rep$legacy_end <- sequence_edta_no_rep$end
      sequence_edta_no_rep$adjustment <- compute_adjustment(sequence_edta_no_rep$start)
      sequence_edta_no_rep$start <- sequence_edta_no_rep$start - sequence_edta_no_rep$adjustment
      sequence_edta_no_rep$end <- sequence_edta_no_rep$end - sequence_edta_no_rep$adjustment
      
      # Peak identification (Memory Optimized)
      te_ranges <- IRanges(start = sequence_edta_no_rep$start, end = sequence_edta_no_rep$end)
      te_ranges <- te_ranges[end(te_ranges) > 0]
      start(te_ranges)[start(te_ranges) < 1] <- 1
      
      if (length(te_ranges) == 0) next
      
      te_cov <- coverage(te_ranges)
      min_c <- min(start(te_ranges))
      max_c <- max(end(te_ranges))
      breaks <- seq(min_c, max_c, length.out = 25)
      
      bin_starts <- floor(breaks[-length(breaks)])
      bin_ends <- floor(breaks[-1]) - 1
      bin_ends[length(bin_ends)] <- floor(breaks[length(breaks)])
      bin_starts[bin_starts < 1] <- 1
      len_cov <- length(te_cov)
      bin_ends[bin_ends > len_cov] <- len_cov
      bin_starts[bin_starts > len_cov] <- len_cov
      
      valid_bins <- bin_starts <= bin_ends
      bin_counts <- numeric(length(bin_starts))
      if(any(valid_bins)) {
        v <- Views(te_cov, start = bin_starts[valid_bins], end = bin_ends[valid_bins])
        bin_counts[valid_bins] <- viewSums(v)
      }
      
      counts <- c(bin_counts[1], bin_counts[1], bin_counts, bin_counts[length(bin_counts)], bin_counts[length(bin_counts)])
      ma_values <- ma(counts)[3 : (length(counts) - 2)]
      mids <- (breaks[-length(breaks)] + breaks[-1]) / 2
      edta_peak <- mids[which.max(ma_values)]
      
      # Distance Distribution (Memory Optimized)
      r_right <- te_ranges[start(te_ranges) >= edta_peak]
      d_right <- IRanges(start = floor(start(r_right) - edta_peak), end = floor(end(r_right) - edta_peak))
      
      r_left <- te_ranges[end(te_ranges) <= edta_peak]
      d_left <- IRanges(start = floor(edta_peak - end(r_left)), end = floor(edta_peak - start(r_left)))
      
      r_cross <- te_ranges[start(te_ranges) < edta_peak & end(te_ranges) > edta_peak]
      d_cross_1 <- IRanges(start = 0, end = floor(edta_peak - start(r_cross)))
      d_cross_2 <- IRanges(start = 0, end = floor(end(r_cross) - edta_peak))
      
      all_dist_ranges <- c(d_right, d_left, d_cross_1, d_cross_2)
      all_dist_ranges <- all_dist_ranges[width(all_dist_ranges) > 0]
      
      if(length(all_dist_ranges) == 0) {
         lm_coef_TE <- list(coefficients = c(0, 0))
      } else {
        dist_cov <- coverage(all_dist_ranges)
        max_dist <- max(end(all_dist_ranges))
        d_breaks <- seq(0, max_dist, length.out = 26)
        
        d_bin_starts <- floor(d_breaks[-length(d_breaks)]) + 1
        d_bin_ends <- floor(d_breaks[-1])
        d_bin_starts[1] <- 1
        len_d_cov <- length(dist_cov)
        d_bin_ends[d_bin_ends > len_d_cov] <- len_d_cov
        
        d_counts <- numeric(length(d_bin_starts))
        valid_d <- d_bin_starts <= d_bin_ends
        if(any(valid_d)) {
          v_d <- Views(dist_cov, start = d_bin_starts[valid_d], end = d_bin_ends[valid_d])
          d_counts[valid_d] <- viewSums(v_d)
        }
        
        d_mids <- (d_breaks[-length(d_breaks)] + d_breaks[-1]) / 2
        d_widths <- diff(d_breaks)
        d_density <- d_counts / d_widths / 2
        
        dist_to_mid_hist <- data.frame(
          dist_to_mid = log10(100 * d_mids / chromosome_no_rep_size),
          counts = d_density
        )
        dist_to_mid_hist <- dist_to_mid_hist[dist_to_mid_hist$counts > 0 & is.finite(dist_to_mid_hist$dist_to_mid),]
        
        if(nrow(dist_to_mid_hist) < 2) {
          lm_coef_TE <- list(coefficients = c(0, 0))
        } else {
          lm_coef_TE <- lm(counts ~ dist_to_mid, dist_to_mid_hist)
        }
      }
    }
    
    if(!no_heli) {
      # calculate gene landscape for -proximity scoring
      sequence_genes <- genes[genes$seqID == chromosomes[j],]
      if(nrow(sequence_genes) == 0) next
      sequence_genes_no_rep <- sequence_genes
      sequence_genes_no_rep$legacy_start <- sequence_genes_no_rep$start
      sequence_genes_no_rep$legacy_end <- sequence_genes_no_rep$end
      sequence_genes_no_rep$adjustment <- compute_adjustment(sequence_genes_no_rep$start)
      sequence_genes_no_rep$start <- sequence_genes_no_rep$start - sequence_genes_no_rep$adjustment
      sequence_genes_no_rep$end <- sequence_genes_no_rep$end - sequence_genes_no_rep$adjustment
      
      # Valley identification (Memory Optimized)
      gene_ranges <- IRanges(start = sequence_genes_no_rep$start, end = sequence_genes_no_rep$end)
      gene_ranges <- gene_ranges[end(gene_ranges) > 0]
      start(gene_ranges)[start(gene_ranges) < 1] <- 1
      
      if (length(gene_ranges) == 0) next
      
      gene_cov <- coverage(gene_ranges)
      min_c <- min(start(gene_ranges))
      max_c <- max(end(gene_ranges))
      breaks <- seq(min_c, max_c, length.out = 25)
      
      bin_starts <- floor(breaks[-length(breaks)])
      bin_ends <- floor(breaks[-1]) - 1
      bin_ends[length(bin_ends)] <- floor(breaks[length(breaks)])
      bin_starts[bin_starts < 1] <- 1
      len_cov <- length(gene_cov)
      bin_ends[bin_ends > len_cov] <- len_cov
      bin_starts[bin_starts > len_cov] <- len_cov
      
      valid_bins <- bin_starts <= bin_ends
      bin_counts <- numeric(length(bin_starts))
      if(any(valid_bins)) {
        v <- Views(gene_cov, start = bin_starts[valid_bins], end = bin_ends[valid_bins])
        bin_counts[valid_bins] <- viewSums(v)
      }
      
      counts <- c(bin_counts[1], bin_counts[1], bin_counts, bin_counts[length(bin_counts)], bin_counts[length(bin_counts)])
      ma_values <- ma(counts)[3 : (length(counts) - 2)]
      mids <- (breaks[-length(breaks)] + breaks[-1]) / 2
      gene_valley <- mids[which.min(ma_values)]
      
      
      # Distance Distribution (Memory Optimized)
      r_right <- gene_ranges[start(gene_ranges) >= gene_valley]
      d_right <- IRanges(start = floor(start(r_right) - gene_valley), end = floor(end(r_right) - gene_valley))
      
      r_left <- gene_ranges[end(gene_ranges) <= gene_valley]
      d_left <- IRanges(start = floor(gene_valley - end(r_left)), end = floor(gene_valley - start(r_left)))
      
      r_cross <- gene_ranges[start(gene_ranges) < gene_valley & end(gene_ranges) > gene_valley]
      d_cross_1 <- IRanges(start = 0, end = floor(gene_valley - start(r_cross)))
      d_cross_2 <- IRanges(start = 0, end = floor(end(r_cross) - gene_valley))
      
      all_dist_ranges <- c(d_right, d_left, d_cross_1, d_cross_2)
      all_dist_ranges <- all_dist_ranges[width(all_dist_ranges) > 0]
      
      if(length(all_dist_ranges) == 0) {
         lm_coef_genes <- list(coefficients = c(0, 0))
      } else {
        dist_cov <- coverage(all_dist_ranges)
        max_dist <- max(end(all_dist_ranges))
        d_breaks <- seq(0, max_dist, length.out = 26)
        
        d_bin_starts <- floor(d_breaks[-length(d_breaks)]) + 1
        d_bin_ends <- floor(d_breaks[-1])
        d_bin_starts[1] <- 1
        len_d_cov <- length(dist_cov)
        d_bin_ends[d_bin_ends > len_d_cov] <- len_d_cov
        
        d_counts <- numeric(length(d_bin_starts))
        valid_d <- d_bin_starts <= d_bin_ends
        if(any(valid_d)) {
          v_d <- Views(dist_cov, start = d_bin_starts[valid_d], end = d_bin_ends[valid_d])
          d_counts[valid_d] <- viewSums(v_d)
        }
        
        d_mids <- (d_breaks[-length(d_breaks)] + d_breaks[-1]) / 2
        d_widths <- diff(d_breaks)
        d_density <- d_counts / d_widths / 2
        
        dist_to_mid_hist <- data.frame(
          dist_to_mid = log10(100 * d_mids / chromosome_no_rep_size),
          counts = d_density
        )
        dist_to_mid_hist <- dist_to_mid_hist[dist_to_mid_hist$counts > 0 & is.finite(dist_to_mid_hist$dist_to_mid),]
        
        if(nrow(dist_to_mid_hist) < 2) {
          lm_coef_genes <- list(coefficients = c(0, 0))
        } else {
          lm_coef_genes <- lm(counts ~ dist_to_mid, dist_to_mid_hist)
        }
      }
      
    }
    
    sequence_repeats$adjustment <- compute_adjustment(sequence_repeats$start)
    sequence_repeats$start_adj <- sequence_repeats$start - sequence_repeats$adjustment
    sequence_repeats$end_adj <- sequence_repeats$end - sequence_repeats$adjustment
    sequence_repeats$start_adj[sequence_repeats$start_adj < 0] = 0
    
    
    # Genes and TEs setup complete, find classes to be scored and initialise the data table using MAX_REPEATS_PER_CHROMOSOME
    
    if(length(top_N_unique)) {
      chr_classes <- data.frame(class = top_N_unique)
      chr_classes$new_class_num_ID <- seq_along(top_N_unique)
      chr_classes$count <- 0
      chr_classes$mean_length <- 0
      chr_classes$total_bp <- 0
      for(k in seq_len(nrow(chr_classes))) {
        chr_classes$count[k] <- nrow(sequence_repeats[sequence_repeats$new_class == chr_classes$class[k], ])
        chr_classes$mean_length[k] <- mean(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
        chr_classes$total_bp[k] <- sum(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
      }
      chr_classes <- chr_classes[order(chr_classes$total_bp, decreasing = TRUE), ]
      chr_classes <- chr_classes[chr_classes$total_bp > 0,]
      if(nrow(chr_classes) > MAX_REPEATS_PER_CHROMOSOME) {
        chr_classes <- chr_classes[1 : MAX_REPEATS_PER_CHROMOSOME, ]
      } 
    } else {
      chr_classes <- data.frame()
    }
    
    
    ### Calculate predictor values ===============================================
    print("predictors")
    if(nrow(chr_classes) == 0) next
    chr_classes$chromosome <- chromosomes[j]
    chr_classes$total_bp_norm_chr <- -1
    chr_classes$total_bp_norm_rep <- -1
    chr_classes$start_sd_norm_chr <- -1
    chr_classes$start_norm_chr_0_50 <- -1
    chr_classes$gaps_count <- -1
    chr_classes$gaps_with_TEs_fraction <- -1
    chr_classes$centre_array_edit <- -1
    chr_classes$centre_array_width_sd <- -1
    chr_classes$centre_chromosome_edit <- -1
    chr_classes$centre_chromosome_width_sd <- -1
    chr_classes$array_sizes_sd_norm_mean_arr_size <- -1
    chr_classes$array_count <- -1
    chr_classes$TE_prox_dist <- -1
    chr_classes$TE_prox_SD <- -1
    chr_classes$TE_lm_coef <- -1
    chr_classes$TE_prox_score <- -1
    chr_classes$gene_prox_dist <- -1
    chr_classes$gene_prox_SD <- -1
    chr_classes$gene_lm_coef <- -1
    chr_classes$gene_prox_score <- -1
    chr_classes$pred_centrophilic_TE <- "LTR/Gypsy"
    chr_classes$t_test_t_val <- -1
    chr_classes$t_test_p_val <- -1
    
    
    
    for(k in seq_len(nrow(chr_classes))) {
      chr_family_repeats <- sequence_repeats[sequence_repeats$new_class == chr_classes$class[k], ]
      chr_family_arrays <- arrays[arrays$new_class_num_ID == chr_classes$new_class_num_ID[k] & arrays$seqID == chromosomes[j], ]
      cat("chromosome no", j, "/", length(chromosomes), "class no", k, "/", nrow(chr_classes), "arrays no:", nrow(chr_family_arrays), "\n")
      # if(nrow(chr_family_arrays) == 0) next
      if(nrow(chr_family_repeats) == 0) next
      # What is the repeat size?
      #   Calculation: Mean repeat length
      chr_classes$mean_length[k] <- chr_classes$mean_length[k] ### PREDICTOR ###
      
      # How big is the family?
      #   Calculation 1: Total bp normalized by chromosome length
      #   Calculation 2: Total bp normalized by all chromosome repeats bp
      chr_classes$total_bp_norm_chr[k] <- chr_classes$total_bp[k] / chromosomes_lengths[j] ### PREDICTOR ###
      chr_classes$total_bp_norm_rep[k] <- chr_classes$total_bp[k] / sum(sequence_repeats$width) ### PREDICTOR ###
      
      # How concentrated is the family? Low score could be holocentric!
      #   Calculation: Start positions SD normalised by chromosome length
      chr_classes$start_sd_norm_chr[k] <- sd(chr_family_repeats$start) / chromosomes_lengths[j] ### PREDICTOR ###
      
      # Where are the repeats?
      #   Calculation: Mean start positions value normalized by chromosome length, presented as distance to the closest chromosome edge (values 0:50)
      chr_classes$start_norm_chr_0_50[k] <- mean(chr_family_repeats$start) / chromosomes_lengths[j] ### PREDICTOR ###
      if(chr_classes$start_norm_chr_0_50[k] > 0.5) chr_classes$start_norm_chr_0_50[k] <- 1 - chr_classes$start_norm_chr_0_50[k] ### PREDICTOR ###
      
      # What are the array interspersed elements (if any)?
      #   Calculation: Find array interspersed sequence (under 20 kbp gaps, more than 200 bp) and count what fraction of them are Tes or TE-derived sequences
      if(!no_edta) {
        if(nrow(chr_family_repeats) > 1) {
          gaps_start <- chr_family_repeats$end[1 : (nrow(chr_family_repeats) - 1)] + 1
          gaps_size <- chr_family_repeats$start[2 : (nrow(chr_family_repeats) - 0)] - chr_family_repeats$end[1 : (nrow(chr_family_repeats) - 1)] - 1
          gaps_start <- gaps_start[gaps_size >= 200 & gaps_size <= 20000]
          gaps_size <- gaps_size[gaps_size >= 200 & gaps_size <= 20000]
          chr_classes$gaps_count[k] <- length(gaps_size) ### PREDICTOR ###
          if(length(gaps_start) != 0) {
            gaps_filled <- rep(FALSE, length(gaps_start))
            for(i_gap in seq_along(gaps_size)) {
              gr3 <- GRanges(chromosomes[j], IRanges(gaps_start[i_gap], (gaps_start[i_gap] + gaps_size[i_gap])))
              
              overlaps <- as.data.frame(findOverlaps(gr3, gr2)) # gr2 are EDTA annotations
              if(nrow(overlaps) != 0) {
                gap_starts <- sequence_edta$start[overlaps$subjectHits]
                gap_ends <- sequence_edta$end[overlaps$subjectHits]
                gap_overlaps <- unlist(lapply(seq_along(gap_starts), function(X) sum( (gap_starts[X] : gap_ends[X]) %in% (gaps_start[i_gap] : (gaps_start[i_gap] + gaps_size[i_gap])) )))
                gap_overlaps_fraction <- gap_overlaps / gaps_size[i_gap]
                if(sum(gap_overlaps_fraction) > 0.25) gaps_filled[i_gap] <- TRUE
              }
              
            }
            chr_classes$gaps_with_TEs_fraction[k] <- sum(gaps_filled) / chr_classes$gaps_count[k] ### PREDICTOR ###
          }
        }
      } 
      
      min_repeats_to_align <- 100
      max_repeats_to_align <- 1000
      desired_fraction_to_align <- 0.2
      
      number_of_repeats_scored <- 0
      cumulative_adist_score <- 0
      mean_centre_width_SD <- NULL
      centre_repeats_count <- NULL
      arrays_sizes <- NULL
      for(array_id in seq_len(nrow(chr_family_arrays))) {
        # What is the repeat divergence within the central parts of the array?
        #   Calculation: central 50% repeats edit distance to THEIR consensus normalized by repeat mean length
        array_repeats <- chr_family_repeats[chr_family_repeats$arrayID == chr_family_arrays$arrayID[array_id], ]
        if(nrow(array_repeats) == 0) {
          arrays_sizes <- c(arrays_sizes, 0)
        } else {
          arrays_sizes <- c(arrays_sizes, sum(array_repeats$width))
        }
        
        if(nrow(array_repeats) < 10) next
        
        repeats_to_align <- round(nrow(array_repeats) * desired_fraction_to_align)
        if(repeats_to_align > max_repeats_to_align) repeats_to_align <- max_repeats_to_align
        if(repeats_to_align < min_repeats_to_align) repeats_to_align <- min_repeats_to_align
        if(repeats_to_align > nrow(array_repeats)) repeats_to_align <- nrow(array_repeats)
        
        repeats_to_align_IDs <- sample(1 : nrow(array_repeats), repeats_to_align)
        
        
        sequences_to_align <- array_repeats$sequence[repeats_to_align_IDs]
        
        # Skip if no valid sequences
        if(length(sequences_to_align) == 0 || all(is.na(sequences_to_align))) next
        
        a <- capture.output({alignment_matrix = suppressWarnings(msa(sequences_to_align, method = "ClustalOmega", type = "dna"))})
        centre_consensus <- consensus_N(alignment_matrix, round(mean(nchar(sequences_to_align))))
        cumulative_adist_score = cumulative_adist_score + sum(adist(centre_consensus, sequences_to_align))
        number_of_repeats_scored <- number_of_repeats_scored + length(sequences_to_align)
        
        # What is the repeat length variation within the central parts of the array?
        #   Calculation: central 50% repeats length SD normalized by repeat mean length
        mean_centre_width_SD <- c(mean_centre_width_SD, sd(nchar(sequences_to_align)))
        centre_repeats_count <- c(centre_repeats_count, length(sequences_to_align))
      }
      # cat("\n")
      if(number_of_repeats_scored != 0) {
        chr_classes$centre_array_edit[k] <- cumulative_adist_score / number_of_repeats_scored ### PREDICTOR ###
        chr_classes$centre_array_width_sd[k] <- sum(mean_centre_width_SD * centre_repeats_count / sum(centre_repeats_count)) ### PREDICTOR ###
      }
      
      if(nrow(chr_family_repeats) > 10) {
        # What is the repeat divergence within the central parts of the chromosome?
        #   Calculation: As above, to identify holocentrics
        
        repeats_to_align <- round(nrow(chr_family_repeats) * desired_fraction_to_align)
        if(repeats_to_align > max_repeats_to_align) repeats_to_align <- max_repeats_to_align
        if(repeats_to_align < min_repeats_to_align) repeats_to_align <- min_repeats_to_align
        if(repeats_to_align > nrow(chr_family_repeats)) repeats_to_align <- nrow(chr_family_repeats)
        
        repeats_to_align_IDs <- sample(1 : nrow(chr_family_repeats), repeats_to_align)
        
        sequences_to_align <- chr_family_repeats$sequence[repeats_to_align_IDs]
        
        # Skip if no valid sequences
        if(length(sequences_to_align) == 0 || all(is.na(sequences_to_align))) next
        
        # cat("Chr ", chromosomes[j], ", class ", chr_classes$class[k],", central chromosome repeats (", nrow(chr_family_repeats), " reps)", sep = "")
        
        # a <- capture.output({alignment_matrix = tolower(as.matrix(msa(sequences_to_align, method = "ClustalOmega", type = "dna")))})
        a <- capture.output({alignment_matrix = suppressWarnings(msa(sequences_to_align, method = "ClustalOmega", type = "dna"))})
        centre_consensus <- consensus_N(alignment_matrix, round(mean(nchar(sequences_to_align))))
        
        chr_classes$centre_chromosome_edit[k] <- sum(adist(centre_consensus, sequences_to_align)) / length(sequences_to_align) ### PREDICTOR ###
        
        # What is the repeat length variation within the central parts of the chromosome?
        #   Calculation: As above, to identify holocentrics
        chr_classes$centre_chromosome_width_sd[k] <- sd(nchar(sequences_to_align)) ### PREDICTOR ###
      }
      
      
      # Can we find the 'main' array?
      #   Calculation: What is the SD of array sizes normalised by mean array size
      chr_classes$array_count[k] <- length(arrays_sizes)
      if(!length(arrays_sizes)) arrays_sizes = 0
      chr_classes$array_sizes_sd_norm_mean_arr_size[k] <- sd(arrays_sizes) / mean(arrays_sizes) ### PREDICTOR ###
      if(length(arrays_sizes) == 1) {
        chr_classes$array_sizes_sd_norm_mean_arr_size[k] = 0
      }
      
      #   What is the TE landscape in the proximity?
      if(!no_edta) {
        chr_classes$TE_prox_dist[k] <- 100 * abs(mean(chr_family_repeats$start_adj) - edta_peak) / chromosome_no_rep_size # the lower the better
        
        chr_classes$TE_prox_SD[k] <- sd(100 * abs(chr_family_repeats$start_adj - edta_peak) / chromosome_no_rep_size) # the lower the better
        
        chr_classes$TE_lm_coef[k] <- lm_coef_TE$coefficients[2] # the higher absolute the better
        
      }
      
      # Score the repeat class based on it’s own mean start position distance to the TE “peak”, score accounting for the correlation calculated before
      
      if(!no_edta) {
        chr_classes$TE_prox_score[k] <- abs(chr_classes$TE_lm_coef[k]) / (chr_classes$TE_prox_dist[k] + chr_classes$TE_prox_SD[k]) * 100
        
      } else {
        chr_classes$TE_prox_score[k] <- 0
      }
      
      # What is the gene landscape in the proximity?
      #   Calculation: Identical to the TE, but this time expect negative correlation
      if(!no_heli) {
        chr_classes$gene_prox_dist[k] <- 100 * abs(mean(chr_family_repeats$start_adj) - gene_valley) / chromosome_no_rep_size # the lower the better
        
        chr_classes$gene_prox_SD[k] <- sd(100 * abs(chr_family_repeats$start_adj - gene_valley) / chromosome_no_rep_size) # the lower the better
        
        chr_classes$gene_lm_coef[k] <- lm_coef_genes$coefficients[2] # the higher the better
        
        chr_classes$gene_prox_score[k] <-  abs(chr_classes$gene_lm_coef[k]) / (chr_classes$gene_prox_dist[k] + chr_classes$gene_prox_SD[k]) * 100
        
      } else {
        chr_classes$gene_prox_score[k] <-  0
      }
      
      
      # What TE families can be identified in the proximity?
      #   Calculation: For each individual TE, calculate distance between it and the closest repeat unit of analysed family, then divide the mean distances of predicted centrophilic TEs by mean distances of other TEs, the lower the ratio, the higher the score
      
      if(!no_edta) {
        sequence_edta$centrophilic_Classification = FALSE
        sequence_edta$centrophilic_Classification[grep(chr_classes$pred_centrophilic_TE[1], sequence_edta$Classification)] = TRUE
        
        sequence_edta$dist_to_closest_rep <- 0
        
        for(seq_edta_id in seq_len(nrow(sequence_edta))) {
          sequence_edta$dist_to_closest_rep[seq_edta_id] <- min(abs(sequence_edta$start[seq_edta_id] - chr_family_repeats$start))
        }
        
        chr_classes$t_test_p_val[k] <- -1
        chr_classes$t_test_t_val[k] <- -1
        if(length(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification]) < 6) next
        if(length(unique(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification])) == 1) next
        chr_classes$t_test_p_val[k] <-  t.test(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification], 
                                               sequence_edta$dist_to_closest_rep[!sequence_edta$centrophilic_Classification])$p.value
        
        chr_classes$t_test_t_val[k] <- t.test(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification], 
                                              sequence_edta$dist_to_closest_rep[!sequence_edta$centrophilic_Classification])$statistic
        
      } 
    }
    # adjust some scores
  
    chr_classes$centre_array_edit <- 100 - (100 * chr_classes$centre_array_edit / chr_classes$mean_length)
    chr_classes$centre_array_width_sd <- 100 - (100 * chr_classes$centre_array_width_sd / chr_classes$mean_length)
    chr_classes$centre_chromosome_edit <- 100 - (100 * chr_classes$centre_chromosome_edit / chr_classes$mean_length)
    chr_classes$centre_chromosome_width_sd <- 100 - (100 * chr_classes$centre_chromosome_width_sd / chr_classes$mean_length)
    
    chr_classes$centre_array_edit[chr_classes$centre_array_edit > 100] = NA
    chr_classes$centre_array_width_sd[chr_classes$centre_array_width_sd > 100] = NA
    chr_classes$centre_chromosome_edit[chr_classes$centre_chromosome_edit > 100] = NA
    chr_classes$centre_chromosome_width_sd[chr_classes$centre_chromosome_width_sd > 100] = NA
    
    chr_classes$centre_array_edit[chr_classes$centre_array_edit < 0] = NA
    chr_classes$centre_array_width_sd[chr_classes$centre_array_width_sd < 0] = NA
    chr_classes$centre_chromosome_edit[chr_classes$centre_chromosome_edit < 0] = NA
    chr_classes$centre_chromosome_width_sd[chr_classes$centre_chromosome_width_sd < 0] = NA
    
    chr_classes$total_bp_norm_chr[chr_classes$total_bp_norm_chr < 0] = NA
    chr_classes$total_bp_norm_rep[chr_classes$total_bp_norm_rep < 0] = NA
    chr_classes$start_sd_norm_chr[chr_classes$start_sd_norm_chr < 0] = NA
    chr_classes$start_norm_chr_0_50[chr_classes$start_norm_chr_0_50 < 0] = NA
    chr_classes$gaps_count[chr_classes$gaps_count < 0] = 0
    chr_classes$gaps_with_TEs_fraction[chr_classes$gaps_with_TEs_fraction < 0] = NA
    
    chr_classes$array_sizes_sd_norm_mean_arr_size[chr_classes$array_sizes_sd_norm_mean_arr_size < 0] = NA
    chr_classes$array_count[chr_classes$array_count < 0] = NA
    chr_classes$TE_prox_dist[chr_classes$TE_prox_dist < 0] = NA
    chr_classes$TE_prox_SD[chr_classes$TE_prox_SD < 0] = NA
    chr_classes$TE_lm_coef[chr_classes$TE_lm_coef == -1] = NA
    chr_classes$TE_prox_score[chr_classes$TE_prox_score < 0] = NA
    chr_classes$gene_prox_dist[chr_classes$gene_prox_dist < 0] = NA
    chr_classes$gene_prox_SD[chr_classes$gene_prox_SD < 0] = NA
    chr_classes$gene_lm_coef[chr_classes$gene_lm_coef == -1] = NA
    chr_classes$gene_prox_score[chr_classes$gene_prox_score < 0] = NA
    chr_classes$t_test_t_val[chr_classes$t_test_t_val == -1] = NA
    chr_classes$t_test_p_val[chr_classes$t_test_p_val == -1] = NA

    all_scores <- rbind(all_scores, chr_classes)
  }
  
  
  
}

# Save output
if(nrow(all_scores) == 0) {
  all_scores <- data.frame(
    class = character(), new_class_num_ID = numeric(), count = numeric(),
    mean_length = numeric(), total_bp = numeric(), chromosome = character(),
    total_bp_norm_chr = numeric(), total_bp_norm_rep = numeric(), 
    start_sd_norm_chr = numeric(), start_norm_chr_0_50 = numeric(),
    gaps_count = numeric(), gaps_with_TEs_fraction = numeric(),
    centre_array_edit = numeric(), centre_array_width_sd = numeric(),
    centre_chromosome_edit = numeric(), centre_chromosome_width_sd = numeric(),
    array_sizes_sd_norm_mean_arr_size = numeric(), array_count = numeric(),
    TE_prox_dist = numeric(), TE_prox_SD = numeric(), TE_lm_coef = numeric(),
    TE_prox_score = numeric(), gene_prox_dist = numeric(), gene_prox_SD = numeric(),
    gene_lm_coef = numeric(), gene_prox_score = numeric(),
    pred_centrophilic_TE = character(), t_test_t_val = numeric(), t_test_p_val = numeric()
  )
}
write.csv(all_scores, file = output_centromeric_scores, row.names = FALSE)


