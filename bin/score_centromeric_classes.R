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


# Load optionals if provided
# print("load edta")
if (!no_edta) edta <- read.csv(tes_filtered_csv)
# print("load genes")
if (!no_heli) genes <- read.csv(genes_filtered_csv)


# Score

{
  chromosomes <- metadata_data$chromosome.name
  chromosomes_lengths <- metadata_data$size
  
  
  ### Decide which classes to score
  
  # Filtering variables
  MIN_BP_TO_SCORE <- 1000
  MIN_REPEAT_COUNT_TO_SCORE <- 3
  MIN_REP_WIDTH_TO_SCORE <- 8
  N_MAX_REP_PER_CHR <- 5
  MAX_REPEATS_PER_CHROMOSOME <- 10
  
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
      sequence_edta = edta[edta$V1 == chromosomes[j], ]
      gr2 <- with(sequence_edta, GRanges(chromosomes[j], IRanges(V4, V5)))
      sequence_edta_no_rep <- sequence_edta
      sequence_edta_no_rep$legacy_V4 <- sequence_edta_no_rep$V4
      sequence_edta_no_rep$legacy_V5 <- sequence_edta_no_rep$V5
      sequence_edta_no_rep$adjustment <- compute_adjustment(sequence_edta_no_rep$V4)
      sequence_edta_no_rep$V4 <- sequence_edta_no_rep$V4 - sequence_edta_no_rep$adjustment
      sequence_edta_no_rep$V5 <- sequence_edta_no_rep$V5 - sequence_edta_no_rep$adjustment
      
      # Peak identification
      TE_coordinates <- list()
      for(edta_id in seq_len(nrow(sequence_edta_no_rep))) {
        TE_coordinates <- append(TE_coordinates, list(sequence_edta_no_rep$V4[edta_id] : sequence_edta_no_rep$V5[edta_id]))
      }
      TE_coordinates <- unlist(TE_coordinates)
      length(TE_coordinates)
      length(unique(TE_coordinates))
      # find arithmetic mean of all scores as the “peak”
      hist_EDTA <- hist(TE_coordinates, breaks = seq(min(TE_coordinates), max(TE_coordinates), length.out = 25), plot = FALSE)
      counts <- c(hist_EDTA$counts[1], hist_EDTA$counts[1], hist_EDTA$counts, hist_EDTA$counts[length(hist_EDTA$counts)], hist_EDTA$counts[length(hist_EDTA$counts)])
      ma_values <- ma(counts)[3 : (length(counts) - 2)]
      edta_peak <- hist_EDTA$mids[which.max(ma_values)]
      
      # calculate TE bp positions density in 2% chromosome length bins
      # calculate window TE density vs distance to peak and test it’s correlation

      distance_to_mid <- abs(edta_peak - TE_coordinates)
      
      dist_to_mid_hist <- data.frame(dist_to_mid = log10(100 * dist_to_mid_hist$mids / chromosome_no_rep_size),
                                     counts = (dist_to_mid_hist$counts) / (dist_to_mid_hist$breaks[2:length(dist_to_mid_hist$breaks)] - dist_to_mid_hist$breaks[1:(length(dist_to_mid_hist$breaks) - 1)]) / 2)
      dist_to_mid_hist <- dist_to_mid_hist[dist_to_mid_hist$counts != -Inf,]
      summary(lm(counts ~ dist_to_mid, dist_to_mid_hist))
      lm_coef_TE <- lm(counts ~ dist_to_mid, dist_to_mid_hist)
      
      
    }
    
    if(!no_heli) {
      # calculate gene landscape for -proximity scoring
      sequence_genes <- genes[genes$V1 == chromosomes[j],]
      sequence_genes_no_rep <- sequence_genes
      sequence_genes_no_rep$legacy_V4 <- sequence_genes_no_rep$V4
      sequence_genes_no_rep$legacy_V5 <- sequence_genes_no_rep$V5
      sequence_genes_no_rep$adjustment <- compute_adjustment(sequence_genes_no_rep$V4)
      sequence_genes_no_rep$V4 <- sequence_genes_no_rep$V4 - sequence_genes_no_rep$adjustment
      sequence_genes_no_rep$V5 <- sequence_genes_no_rep$V5 - sequence_genes_no_rep$adjustment
      
      
      # Valley identification
      gene_coordinates <- list()
      for(edta_id in seq_len(nrow(sequence_genes_no_rep))) {
        gene_coordinates <- append(gene_coordinates, list(sequence_genes_no_rep$V4[edta_id] : sequence_genes_no_rep$V5[edta_id]))
      }
      gene_coordinates <- unlist(gene_coordinates)
      length(gene_coordinates)
      length(unique(gene_coordinates))
      
      hist_gene <- hist(gene_coordinates, breaks = seq(min(gene_coordinates), max(gene_coordinates), length.out = 25), plot = FALSE)
      counts <- c(hist_gene$counts[1], hist_gene$counts[1], hist_gene$counts, hist_gene$counts[length(hist_gene$counts)], hist_gene$counts[length(hist_gene$counts)])
      ma_values <- ma(counts)[3 : (length(counts) - 2)]
      gene_valley <- hist_gene$mids[which.min(ma_values)]
      
      
      distance_to_mid <- abs(edta_peak - gene_coordinates)
      
      dist_to_mid_hist <- data.frame(dist_to_mid = log10(100 * dist_to_mid_hist$mids / chromosome_no_rep_size),
                                     counts = (dist_to_mid_hist$counts) / (dist_to_mid_hist$breaks[2:length(dist_to_mid_hist$breaks)] - dist_to_mid_hist$breaks[1:(length(dist_to_mid_hist$breaks) - 1)]) / 2)
      dist_to_mid_hist <- dist_to_mid_hist[dist_to_mid_hist$counts != -Inf,]
      summary(lm(counts ~ dist_to_mid, dist_to_mid_hist))
      lm_coef_genes <- lm(counts ~ dist_to_mid, dist_to_mid_hist)
      
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
      cat(k, "/", nrow(chr_classes), "\n")
      chr_family_repeats <- sequence_repeats[sequence_repeats$new_class == chr_classes$class[k], ]
      chr_family_arrays <- arrays[arrays$new_class_num_ID == chr_classes$new_class_num_ID[k] & arrays$seqID == chromosomes[j], ]
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
                gap_starts <- sequence_edta$V4[overlaps$subjectHits]
                gap_ends <- sequence_edta$V5[overlaps$subjectHits]
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
      cat(nrow(chr_family_arrays), "arrays:")
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
        central_repeats <- array_repeats[ceiling(nrow(array_repeats) / 4) : floor(3 * nrow(array_repeats) / 4), ]
        
        repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), round(nrow(central_repeats) * desired_fraction_to_align))
        if(length(repeats_to_align_IDs) > max_repeats_to_align) {
          repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), max_repeats_to_align)
        } else if(length(repeats_to_align_IDs) < min_repeats_to_align) {
          if(nrow(central_repeats) > min_repeats_to_align) {
            repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), min_repeats_to_align)
          } else {
            repeats_to_align_IDs <- seq_len(nrow(central_repeats))
          }
        }
        
        sequences_to_align <- central_repeats$sequence[repeats_to_align_IDs]
        cat(array_id, " (", nrow(array_repeats), " reps), ", sep = "")
        if(length(sequences_to_align) < 2) {
          cat("\n\n\n\n")
          stop(paste0(" did not find repeats in one of the classes: ", classes$class[j], ", investigate"))
        }
        a <- capture.output({alignment_matrix = suppressWarnings(msa(sequences_to_align, method = "ClustalOmega", type = "dna"))})
        centre_consensus <- consensus_N(alignment_matrix, round(mean(nchar(sequences_to_align))))
        cumulative_adist_score = cumulative_adist_score + sum(adist(centre_consensus, sequences_to_align))
        number_of_repeats_scored <- number_of_repeats_scored + length(sequences_to_align)
        
        # What is the repeat length variation within the central parts of the array?
        #   Calculation: central 50% repeats length SD normalized by repeat mean length
        mean_centre_width_SD <- c(mean_centre_width_SD, sd(nchar(sequences_to_align)))
        centre_repeats_count <- c(centre_repeats_count, length(sequences_to_align))
      }
      cat("\n")
      if(number_of_repeats_scored != 0) {
        chr_classes$centre_array_edit[k] <- cumulative_adist_score / number_of_repeats_scored ### PREDICTOR ###
        chr_classes$centre_array_width_sd[k] <- sum(mean_centre_width_SD * centre_repeats_count / sum(centre_repeats_count)) ### PREDICTOR ###
      }
      
      if(nrow(chr_family_repeats) > 10) {
        # What is the repeat divergence within the central parts of the chromosome?
        #   Calculation: As above, to identify holocentrics
        central_repeats <- chr_family_repeats[ceiling(nrow(chr_family_repeats) / 4) : floor(3 * nrow(chr_family_repeats) / 4), ]
        repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), round(nrow(central_repeats) * desired_fraction_to_align))
        if(length(repeats_to_align_IDs) > max_repeats_to_align) {
          repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), max_repeats_to_align)
        } else if(length(repeats_to_align_IDs) < min_repeats_to_align) {
          if(nrow(central_repeats) > min_repeats_to_align) {
            repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), min_repeats_to_align)
          } else {
            repeats_to_align_IDs <- seq_len(nrow(central_repeats))
          }
        }
        sequences_to_align <- central_repeats$sequence[repeats_to_align_IDs]
        cat("Chr ", chromosomes[j], ", class ", chr_classes$class[k],", central chromosome repeats (", nrow(chr_family_repeats), " reps)", sep = "")
        
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
          sequence_edta$dist_to_closest_rep[seq_edta_id] <- min(abs(sequence_edta$V4[seq_edta_id] - chr_family_repeats$start))
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
  
  
}

# Save output
write.csv(chr_classes, file = output_centromeric_scores, row.names = FALSE)  # Assuming scores_data is created in your code


























