#!/usr/bin/env Rscript

# CAP.R
# Description: Generate final CAP outputs with optional TEs and genes

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9 || length(args) > 11) {
  stop("Usage: CAP.R <predictions.csv> <repeats_reclassed.csv> <arrays_reclassed.csv> <genome_classes.csv> <metadata.csv> <assembly_name> <GC> <CTW> [<TEs_filtered.csv>] [<genes_filtered.csv>] [<scores.csv>]")
}
no_edta <- FALSE
no_heli <- FALSE

predictions_csv <- args[1] # ML predictions to match scores_csv
repeats_reclassed_csv <- args[2]
arrays_reclassed_csv <- args[3]
genome_classes_csv <- args[4] # individual classes
metadata_csv <- args[5]
assembly_name <- basename(args[6])
gc_csv <- basename(args[7])
ctw_csv <- basename(args[8])
tes_filtered_csv <- if (args[9] != "NO_FILE") args[9] else no_edta <- TRUE
genes_filtered_csv <- if (args[10] != "NO_FILE") args[10] else no_heli <- TRUE
scores_csv <- basename(args[11]) # scores per class per chromosome

# kmer_data_files <- list.files(path = "/home/pwlodzimierz/ToL/kmer_profiles", full.names = TRUE)


# Load libraries
suppressMessages({
  library(seqinr)
  library(stringr)
  library(Biostrings)
  library(GenomicRanges)
  library(scales)
  # library(BCT)
})

# Load additional functions
source(file.path(Sys.getenv("WORKFLOW_DIR"), "bin", "auxfuns.R"))

# Load required data
print("load main")
predictions_data <- read.csv(predictions_csv)
repeats_data <- read.csv(repeats_reclassed_csv)
arrays_data <- read.csv(arrays_reclassed_csv)
genome_classes_data <- read.csv(genome_classes_csv)
metadata_data <- read.csv(metadata_csv)
gc_data <- read.csv(gc_csv)
ctw_data <- read.csv(ctw_csv)
scores_data <- read.csv(scores_csv)

# Ensure critical data is not empty
if (nrow(metadata_data) == 0) stop("metadata_csv is empty")
if (nrow(genome_classes_data) == 0) genome_classes_data <- data.frame(class = character(), importance = numeric())

# Load optionals if provided
print("load extra")
if (!no_edta) {
  if (file.size(tes_filtered_csv) > 1) {
    tes_data <- read.csv(tes_filtered_csv)
    if (nrow(tes_data) == 0) {
      warning("TEs file is empty, setting no_edta=TRUE")
      no_edta <- TRUE
    }
  } else {
    warning("TEs file is empty (0 bytes), setting no_edta=TRUE")
    no_edta <- TRUE
  }
}
if (!no_heli) {
  if (file.size(genes_filtered_csv) > 1) {
    genes_data <- read.csv(genes_filtered_csv)
    if (nrow(genes_data) == 0) {
      warning("Genes file is empty, setting no_heli=TRUE")
      no_heli <- TRUE
    }
  } else {
    warning("Genes file is empty (0 bytes), setting no_heli=TRUE")
    no_heli <- TRUE
  }
}



repeats  <- repeats_data
arrays   <- arrays_data
classes  <- genome_classes_data
chr_info <- metadata_data
if (nrow(repeats) == 0) repeats <- data.frame(seqID = character(), start = numeric(), width = numeric(), new_class = character())
if (nrow(classes)) classes$num_ID <- seq_len(nrow(classes))

if (!no_heli) {
  genes <- genes_data
  genes <- genes[genes$type == "CDS", ]
}

if (!no_edta) {
  edta <- tes_data
  if ("reassigned" %in% colnames(edta)) {
    # names(edta) <- c("", "V1","V2","V3","V4","V5","V6","V7","V8","ID","Name",
    #                  "Classification","Sequence_ontology","Identity","Method",
    #                  "TSD","TIR","motif","tsd","oldV3","overlapping_bp",
    #                  "width","overlapping_percentage","reassigned")
    edta$start <- edta$start + 1
    edta$end <- edta$end + 1
  }
  edta$start[edta$start == 0] <- 1
}


# ------------------------------------------------------------------ #

####
edta_classes <- list(
  # class I (retrotransposons)
  ## LTR retrotransposons
  c("Gypsy_LTR_retrotransposon"),
  c("Copia_LTR_retrotransposon"),
  c("Bel_Pao_LTR_retrotransposon", "TRIM_LTR_retrotransposon", "Caulimoviridae", "pararetrovirus"),
  c("Retrovirus", "LTR_retrotransposon", "long_terminal_repeat", 
    "Retrovirus_LTR_retrotransposon", "Endogenous_Retrovirus_LTR_retrotransposon"),
  ## Non-LTR retrotransposons
  c("LINE_element", "L1_LINE_retrotransposon", "RTE_LINE_retrotransposon", "Jockey_LINE_retrotransposon", "L2_LINE_retrotransposon", "R1_LINE_retrotransposon", 
    "CR1_LINE_retrotransposon", "I_LINE_retrotransposon", "Rex_LINE_retrotransposon", "Proto2_LINE_retrotransposon", "R2_LINE_retrotransposon", "Tad1_LINE_retrotransposon", "R4_LINE_retrotransposon"),
  c("SINE_element", "tRNA_SINE_retrotransposon", "5S_SINE_retrotransposon", "MIR_SINE_retrotransposon", "Alu_SINE_retrotransposon"),
  c("Penelope_retrotransposon", "DIRS_YR_retrotransposon"),
  c("non_LTR_retrotransposon"),
  # class II (DNA transposons)
  ## TIRs
  c("Kolobok_TIR_transposon" , "Ginger_TIR_transposon", "Academ_TIR_transposon", 
    "Novosib_TIR_transposon", "Sola_TIR_transposon", "Merlin_TIR_transposon", 
    "IS3EU_TIR_transposon", "PiggyBac_TIR_transposon", "hAT_TIR_transposon", 
    "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon", "Dada_TIR_transposon", 
    "CACTA_TIR_transposon", "Zisupton_TIR_transposon", "PIF_Harbinger_TIR_transposon", 
    "P_TIR_transposon", "Transib_TIR_transposon", "PILE_TIR_transposon",
    "POLE_TIR_transposon", "Zator_TIR_transposon"),
  ## other class II
  c("DNA_transposon", "helitron", "MITE"),
  c("Maverick_Polinton", "polinton"),
  # other, recombinase element based
  c("Tyrosine_Recombinase_Elements", "Crypton_Tyrosine_Recombinase", "Crypton_YR_transposon", "Ngaro_YR_retrotransposon"),
  # others
  c("TE", "TE_unclass"),
  # likely not TEs, remove for plotting?
  c("repeat_region", "SUPER", "Sequence_Ontology", "rRNA_gene", "target_site_duplication", "chr", "snRNA", "repeat_fragment"))
edta_classes_colours <-  c(
  "#E31A1C",  # red
  "#D55E00",  # reddish-orange
  "#F5793A",  # bright orange
  # "#FF7F00",  # orange
  # "#FDBF6F",  # peach
  "#F0E442",  # yellow
  # "#6A3D9A",  # dark purple
  "#A95AA1",  # purple
  "#CC79A7",  # pink
  "#DDA0DD",  # light pinkish purple (plum)
  "#CAB2D6",  # lavender
  "#1F78B4",  # blue
  "#B2DF8A",  # light green
  # "#33A02C",  # green
  "#009E73",  # teal green
  # "#56B4E9",  # light blue
  "#7F7F7F",  # grey
  "#000000",   # black
  "#000000"   # black
)

# ------------------------------------------------------------------ #
#  CHROMOSOME FILTERING (top 30 longest)
# ------------------------------------------------------------------ #
# message("Removed from plotting chromosomes shorter than 200 kbp: ", chr_info$chromosome.name[chr_info$size < 200000])
# chr_info <- chr_info[chr_info$size > 200000,] # TODO: remove this, user should provide fasta with only relevant chromosomes or add parameter
chromosomes     <- chr_info$chromosome.name[chr_info$is.chr == 1]
chromosomes_len <- chr_info$size[chr_info$is.chr == 1]

if (length(chromosomes) == 0) stop("No chromosomes found in metadata")
chromosomes_sets <- vector("list", (length(chromosomes)%/%20+1))
chromosomes_len_sets <- vector("list", (length(chromosomes)%/%20+1))

for(j in seq_along(chromosomes)) {
  chromosomes_sets[[(j%/%20)+1]] <- c(chromosomes_sets[[(j%/%20)+1]], chromosomes[j])
  chromosomes_len_sets[[(j%/%20)+1]] <- c(chromosomes_len_sets[[(j%/%20)+1]], chromosomes_len[j])
}

# ------------------------------------------------------------------ #
#  SCORES & CLASS FILTERING
# ------------------------------------------------------------------ #
# scores_file <- file.path(getwd(), paste0("genome_classes_", assembly_name, ".csv"))
scores <- if (nrow(scores_data) > 0 && nrow(predictions_data) > 0) cbind(scores_data, predictions_data) else data.frame()

if (nrow(scores)) {
  # scores$ed_perc <- 100 * scores$centre_array_edit / scores$mean_length
  # scores$width_sd_perc <- 100 * scores$centre_array_width_sd / scores$mean_length
  # scores$ed_perc[scores$ed_perc < 0] = NA
  # scores$width_sd_perc[scores$width_sd_perc < 0] = NA
  # scores$ed_perc = 100 - scores$ed_perc
  # scores$width_sd_perc = 100 - scores$width_sd_perc
  scores$ed_perc <- scores$centre_chromosome_edit
  scores$width_sd_perc <- scores$centre_chromosome_width_sd
}

# ------------------------------------------------------------------ #
#  SELECT CLASSES TO PLOT (top-scoring, >=5 kbp, <=20)
# ------------------------------------------------------------------ #

classes_to_plot <- NULL
if(nrow(scores)) {
  scores$probability_centromeric <- scores$probability_centromeric * 100
  scores <- scores[order(scores$probability_centromeric, decreasing = TRUE), ]
  # scores <- scores[scores$probability_centromeric >= 0.01, ] # optional filtering
  for(j in unique(scores$chromosome)) {
    chr_scores <- scores[scores$chromosome == j, ]
    if(nrow(chr_scores) == 0) next
    if(nrow(chr_scores) > 4) chr_scores <- chr_scores[1:4, ]
    classes_to_plot <- c(classes_to_plot, chr_scores$class[seq_len(min(4, nrow(chr_scores)))])
  }
} 
classes_to_plot <- unique(classes_to_plot)
classes_to_score <- classes$importance[classes$class %in% classes_to_plot]
classes_to_plot <- classes_to_plot[order(classes_to_score, decreasing = T)]
cat("Classes to plot:", classes_to_plot, "\n")


# ------------------------------------------------------------------ #
#  INITIALIZE DENSITY DATA STRUCTURE
# ------------------------------------------------------------------ #
# This will store pre-computed densities for replotting instead of raw data
chromosome_densities <- list()

# ------------------------------------------------------------------ #
#  PLOT NAME
# ------------------------------------------------------------------ #
suffix <- paste0(
  ifelse(no_edta,  "_noedta", "_edta"),
  ifelse(no_heli,  "_nogene", "_gene")
)


# ------------------------------------------------------------------ #
#  PLOTS LOOP SETUP
# ------------------------------------------------------------------ #
cat("Plotting. plots to complete:", length(chromosomes_sets), "\n")

for(k in seq_along(chromosomes_sets)) {
  
  plot_name <- file.path(paste0(assembly_name, "_CAP_plot_", k, "_", suffix, ".png"))
  
  chromosomes <- chromosomes_sets[[k]]
  chromosomes_len <- chromosomes_len_sets[[k]]

  # get kmer data if available
  add_kmers <- FALSE
  kmers_data <- data.frame()
  #if(sum(grepl(assembly_name, kmer_data_files))) {
  #  add_kmers <- TRUE
  #  kmers_data <- read.table(file = kmer_data_files[grep(assembly_name, kmer_data_files)], header = T, sep = "")
  #}
  
  # ------------------------------------------------------------------ #
  #  PLOT SETUP
  # ------------------------------------------------------------------ #
  bin_rep <- 10000
  bin_gc <- 2000 # needs to be the same as in GC.R
  bin_edta <- 100000 # is also used in CTW.R
  bin_gene <- 50
  x_tick  <- 1e6
  palette <- rep(c("#88CCEE","#CC6677","#DDCC77","#117733","#332288",
                   "#AA4499","#44AA99","#999933","#882255","#661100",
                   "#6699CC","#888888","#55AA55","#EE8866","#771155",
                   "#99DDFF","#FFAABB","#4477AA","#D4AABB","#33BBEE"), 100)
  
  png(plot_name, width = 3000, 
      height = (700 + 1200*length(chromosomes)), pointsize = 32)
  
  layout(matrix(1:(4*length(chromosomes)+1), nrow = (1 + 4*length(chromosomes)), ncol = 1), 
         heights = c(700, rep(c(200,400,300,300),length(chromosomes))))
  
  par(mar = c(0.5,2,0.2,0.2), mgp = c(0.5, 0.5, 0), oma = c(2, 3, 3, 4))
  cex_factor <- 1
  # ------------------------------------------------------------------ #
  #  TITLE
  # ------------------------------------------------------------------ #
  plot(NA, xlim = c(1,100), ylim = c(1,100), axes = FALSE, xlab = "", ylab = "")
  text(50, 95, assembly_name, pos = 1, cex = cex_factor*5)
  text(1, 70, sprintf("Total size: %s Mbp", formatC(sum(chromosomes_len)/1e6, format="f", digits=3, big.mark=" ")), pos=4, cex = cex_factor*3)
  text(1,55, sprintf("Chromosomes: %d", length(chromosomes)), pos=4, cex = cex_factor*3)
  text(1,40, sprintf("Transposable elements legend:"), pos=4, cex = cex_factor*2)
  
  text(x = c(1,15,30,45,75), y = 32, 
       labels =  c("class I LTR: ",
                   "Gypsy",
                   "Copia",
                   "Bel Pao;TRIM;Caulimoviridae",
                   "unspecified"), 
       col = c("black", edta_classes_colours[1:4]),
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15,30,45,75), y = 26, 
       labels = c("class I non-LTR: ", 
                  "LINE",
                  "SINE",
                  "Penelope;DIRS YR",
                  "unspecified"),
       col = c("black", edta_classes_colours[5:8]),
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15), y = 20, 
       labels = c("class II TIRs: ", "Kolobok;Ginger;Academ;Novosib;Sola;Merlin;IS3EU;PiggyBac;hAT;Mutator;Tc1 Mariner;Dada;CACTA;Zisupton;PIF Harbinger"),
       col = c("black", edta_classes_colours[9]), 
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15,45), y = 14, 
       labels = c("class II others: ",  
                  "DNA_transposon;helitron;MITE",
                  "Maverick Polinton"),
       col = c("black", edta_classes_colours[10:11]), 
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15,30), y = 8, 
       labels = c("other: ", "Tyrosine Recombinase",
                  "unspecified"),
       col = c("black", edta_classes_colours[12:13]), 
       pos = 4, cex = cex_factor* 1.5)
  
  # ------------------------------------------------------------------ #
  #  PER-CHROMOSOME LOOP
  # ------------------------------------------------------------------ #
  
  for (j in seq_along(chromosomes)) {
    chr <- chromosomes[j]
    len <- chromosomes_len[j]
    rep_chr <- subset(repeats, seqID == chr)
    if (nrow(rep_chr)) {
      rep_chr <- rep_chr[rep_chr$start >= 1,]
      rep_chr <- rep_chr[rep_chr$start <= len,]
    }
    edt_chr <- if (!no_edta && nrow(edta) > 0) subset(edta, seqID == chr) else data.frame()
    gen_chr <- if (!no_heli && nrow(genes) > 0) subset(genes, seqID == chr) else data.frame()
    
    # Initialize density storage for this chromosome
    chr_density_data <- list(
      chromosome = chr,
      length = len,
      rep_coverage_windows = numeric(),
      rep_coverage_values = numeric(),
      edta_coverage_windows = numeric(),
      edta_coverage_values = numeric(),
      family_coverage = list(),
      edta_family_coverage = list(),
      te_repeat_density = numeric(),
      te_repeat_density_mids = numeric(),
      gene_density = numeric(),
      gene_density_mids = numeric()
    )
    
    print(paste0("Genome ", assembly_name, " | Chromosome ", j, "/", length(chromosomes)))
    
    plot(NA,NA, xlim = c(1,100), ylim = c(1,300), xlab = "", ylab = "", axes = F)
    text(x = 1, y = 15, labels = assembly_name, pos = 4, cex = cex_factor* 1.4)
    text(x = 15, y = 15, labels = chromosomes[j], pos = 4, cex = cex_factor* 1.2)
    text(x = 30, y = 15, 
         labels = paste0(formatC(chromosomes_len[j]/1000000, format = "f", big.mark = " ", digits = 3), " Mbp"), 
         pos = 4, cex = cex_factor* 1.2)
    abline(h = 50, lwd = 4)
    abline(h = 60, lwd = 4)
    
    
    # === TABLE ===
    if (nrow(scores) > 0) {
      sc <- subset(scores, chromosome == chr)[, c("class","count","mean_length","total_bp",
                                                  "ed_perc","width_sd_perc","probability_centromeric")]
      sc <- sc[sc$class %in% classes_to_plot,]
      sc[, 3:6] <- round(sc[, 3:6], 2)
      sc[, 7] <- round(sc[, 7], 2)
      sc$colours <- palette[match(sc$class, classes_to_plot)]
      create_table(sc[,1:7], c("Class","Repeats no","Mean width, bp","Total bp",
                               "Sequence similarity %","Width similarity %","Cen probability"),
                   colours = sc$colours,
                   font_size = cex_factor* 1.2)
    } else {
      plot.new()
    }
    
    # === PLOT A: Repeats + GC + families ===
    win_rep <- genomic.bins.starts(1, len, bin.size = bin_rep)
    rep_cov <- if (nrow(rep_chr)) {
      calculate.repeats.percentage.in.windows(win_rep, rep_chr$start, rep_chr$width, len, unique_coords = TRUE)
    } else rep(0, length(win_rep))
    rep_cov[rep_cov == 0] <- NA; rep_cov[rep_cov > 100] <- 100
    
    # Save repeat coverage density
    chr_density_data$rep_coverage_windows <- win_rep
    chr_density_data$rep_coverage_values <- rep_cov
    
    plot(NA, NA, col="#CCCCCC", ylim=c(0,100), xaxt="n", yaxt="n", xlab="", ylab="", xlim = c(0,len))
    mtext("REP% per 10 Kbp", side = 2, line = 0, col = "grey", cex = cex_factor* 0.5, at = 50, adj = 0.5)
    axis(side = 2, labels = c("0","100"), at = c(0,100))

    
    for(window_plot in seq_along(win_rep)) {
      rect(xleft = win_rep[window_plot], ybottom = 0, xright = win_rep[window_plot] + bin_rep, ytop = rep_cov[window_plot], col = "grey", border = NA)
    }
    
    # GC (TODO only if requested)
    gc_chs_data <- gc_data[gc_data$chromosome == chr,]
    if (nrow(gc_chs_data) > 0) {
      gc_mids <- gc_chs_data$bin_mid
      gc_vals <- gc_chs_data$bin_value
      lines(gc_mids, gc_vals, type = "l", lwd = 2, cex = cex_factor * 1)
      mtext("GC% per  2 Kbp", side = 2, line = 2, col = "black", cex = cex_factor* 0.5, at = 50, adj = 0.5)
    }
    
    # CTW (TODO only if requested)
    ctw_chs_data <- ctw_data[ctw_data$chromosome == chr,]
    if (nrow(ctw_chs_data) > 0) {
      ctw_mids <- ctw_chs_data$bin_mid
      ctw_vals <- ctw_chs_data$bin_value
      lines(x = ctw_mids, ctw_vals, type = "l", col = "#FFA500", lwd = 2, cex = cex_factor * 1)
      mtext("CTW per  2 Kbp D10", side = 2, line = 3, col = "#FFA500", cex = cex_factor* 0.5, at = 50, adj = 0.5)
    }
    
    
    # Families
    if (length(classes_to_plot)) {
      for (m in seq_len(length(classes_to_plot))) {
        fam <- subset(rep_chr, new_class == classes_to_plot[m])
        if (nrow(fam) == 0) next
        fam <- fam[fam$start > 0,]
        fam <- fam[fam$end <= len,]
        if(nrow(fam) == 0) next
        cov <- calculate.repeats.percentage.in.windows(win_rep, fam$start, fam$width, len, unique_coords = TRUE)
        cov[cov == 0] <- NA; cov[cov > 100] <- 100
        # Save family coverage
        chr_density_data$family_coverage[[classes_to_plot[m]]] <- cov
        lines(win_rep, cov, col = palette[m], pch=16, type="o", lwd = 2, cex = cex_factor * 1)
      }
      mtext("SIG REP% per 10 Kbp", side = 2, line = 1, col = "#88CCEE", cex = cex_factor* 0.5, at = 50, adj = 0.5)
    }
    axis(1, at = seq(0, len, by = x_tick), labels = FALSE)
    # text(1, bin_rep, chr, pos = 4)
    
    # === PLOT B: EDTA + TE/repeat peak + genes + HiC ===
    win_edta <- genomic.bins.starts(1, len, bin.size = bin_edta)
    edt_chr <- edt_chr[edt_chr$start > 0,]
    edt_chr <- edt_chr[edt_chr$end <= len,]
    edt_cov <- if (!no_edta && nrow(edt_chr) > 0) {
      calculate.repeats.percentage.in.windows(win_edta, edt_chr$start, edt_chr$width, len, unique_coords = TRUE)
    } else rep(0, length(win_edta))
    edt_cov[edt_cov == 0] <- NA; edt_cov[edt_cov > 100] <- 100
    
    # Save EDTA coverage density
    chr_density_data$edta_coverage_windows <- win_edta + bin_edta/2
    chr_density_data$edta_coverage_values <- edt_cov
    
    plot(NA, NA, type="h", col="#CCCCCC", ylim=c(0,100), xaxt="n", yaxt="n", xlab="", ylab="", xlim = c(0,len))
    mtext("EDTA% per 100 Kbp", side = 2, line = 0, col = "grey", cex = cex_factor* 0.5, at = 50, adj = 0.5)
    axis(2, col="grey", at = c(0,100), labels = c("0", "100"))
    
    # add kmer profile background if available
    if(add_kmers) {
      kmers_data_chr <- kmers_data[kmers_data$seqid_or_ == chr,]
      
      kmers_data_chr$bin_start <- kmers_data_chr$dist_start_ * kmers_data_chr$binwidth_
      kmers_data_chr$bin_mid <- kmers_data_chr$bin_start + kmers_data_chr$binwidth_/2
      kmers_data_chr$bin_end <- kmers_data_chr$bin_start + kmers_data_chr$binwidth_
      getout = FALSE
      if(!(1 %in% kmers_data_chr$cutoff_)) {
        message(paste0("Kmer data for ", chr, " does not have cutoff 1, skipping kmer plotting for this chromosome"))
        getout = TRUE
      }
      if(!(25 %in% kmers_data_chr$cutoff_)) {
        message(paste0("Kmer data for ", chr, " does not have cutoff 25, skipping kmer plotting for this chromosome"))
        getout = TRUE
      }
      if(!(50 %in% kmers_data_chr$cutoff_)) {
        message(paste0("Kmer data for ", chr, " does not have cutoff 50, skipping kmer plotting for this chromosome"))
        getout = TRUE
      }
      if(!(75 %in% kmers_data_chr$cutoff_)) {
        message(paste0("Kmer data for ", chr, " does not have cutoff 75, skipping kmer plotting for this chromosome"))
        getout = TRUE
      }
      if(!(90 %in% kmers_data_chr$cutoff_)) {
        message(paste0("Kmer data for ", chr, " does not have cutoff 90, skipping kmer plotting for this chromosome"))
        getout = TRUE
      }
      if(!getout) {
        kmers_data_1 <- kmers_data_chr[kmers_data_chr$cutoff_ == 1,]
        kmers_data_25 <- kmers_data_chr[kmers_data_chr$cutoff_ == 25,]
        kmers_data_50 <- kmers_data_chr[kmers_data_chr$cutoff_ == 50,]
        kmers_data_75 <- kmers_data_chr[kmers_data_chr$cutoff_ == 75,]
        kmers_data_90 <- kmers_data_chr[kmers_data_chr$cutoff_ == 90,]
        
        kmers_data_averaged <- unlist(lapply(
          seq_len(nrow(kmers_data_1)),
          function(X) {
            mean(c(kmers_data_1$occupancy_freq_[X],
                  kmers_data_25$occupancy_freq_[X],
                  kmers_data_50$occupancy_freq_[X],
                  kmers_data_75$occupancy_freq_[X],
                  kmers_data_90$occupancy_freq_[X]))
          }))
        kmers_data_averaged_100kb <- (kmers_data_averaged[1 : (length(kmers_data_averaged) - 1)] + kmers_data_averaged[2 : length(kmers_data_averaged)]) / 2
        if (length(kmers_data_averaged_100kb) > 0) {
          kmers_data_averaged_100kb <- kmers_data_averaged_100kb[seq(1, length(kmers_data_averaged_100kb), by = 2)]
        }
        kmers_data_averaged_100kb <- c(kmers_data_averaged_100kb, kmers_data_averaged[length(kmers_data_averaged)])
        
        if(length(kmers_data_averaged_100kb) != 0) {
          color_func <- colorRamp(c("white", "yellow", "orange", "red"))
          rgb_vals <- color_func(kmers_data_averaged_100kb)
          if(length(rgb_vals) == 0) next
          # Replace NA values in rgb_vals with white (255, 255, 255)
          rgb_vals[is.na(rgb_vals)] <- 255
          colors <- rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue = 255, alpha = 150)
          for(window_plot in seq_along(win_edta)) {
            rect(xleft = win_edta[window_plot] - 50000, ybottom = 0, xright = win_edta[window_plot] + 50000, ytop = 100, col = colors[window_plot], border = NA)
          }
        }
      }
      
      
    }
    for(window_plot in seq_along(win_edta)) {
      rect(xleft = win_edta[window_plot] - 50000, ybottom = 0, xright = win_edta[window_plot] + 50000, ytop = edt_cov[window_plot], col = "#CCCCCC", border = NA)
    }

    # EDTA classes
    if (!no_edta && nrow(edt_chr)) {
      for (m in rev(seq_along(edta_classes))) {
        cls <- edt_chr[edt_chr$type %in% edta_classes[[m]], ]
        cls <- cls[cls$start > 0,]
        cls <- cls[cls$end <= len,]
        if (!nrow(cls)) next
        cov <- calculate.repeats.percentage.in.windows(win_edta, cls$start, cls$width, len, unique_coords = TRUE)
        cov[cov == 0] <- NA; cov[cov > 100] <- 100
        # Save EDTA family coverage
        chr_density_data$edta_family_coverage[[m]] <- list(class_index = m, values = cov)
        lines(win_edta + (bin_edta/2), cov, col = edta_classes_colours[m], pch=16, type="o", lwd = 2, cex = cex_factor * 1)
      }
      mtext("FAM EDTA% per 100 Kbp", side = 2, line = 1, col = "red", cex = cex_factor* 0.5, at = 50, adj = 0.5)
    }
    
    # TE+repeat peak
    te_coords <- c()
    if (!no_edta && nrow(edt_chr)) te_coords <- c(te_coords, unlist(mapply(`:`, edt_chr$start, edt_chr$end)))
    if (nrow(rep_chr))          te_coords <- c(te_coords, unlist(mapply(`:`, rep_chr$start, rep_chr$end)))
    if (length(te_coords) > 0) {
      te_coords <- te_coords[te_coords <= len & te_coords >= 1]
      te_coords <- unique(te_coords)
      if (length(te_coords) > 0) {
        te_hist <- hist(te_coords, breaks = seq(0, len, length.out = bin_gene), plot = FALSE)
        # Only apply moving average if there are enough data points (need at least 5 for padding)
        if (length(te_hist$counts) >= 5) {
          te_ma   <- ma(c(te_hist$counts[1], te_hist$counts[1], te_hist$counts,
                          te_hist$counts[length(te_hist$counts)], te_hist$counts[length(te_hist$counts)]))[3:(length(te_hist$counts)+2)]
          # Save TE+repeat density
          chr_density_data$te_repeat_density <- te_ma
          chr_density_data$te_repeat_density_mids <- te_hist$mids
          par(new = TRUE)
          max_te_ma <- max(te_ma, na.rm = TRUE)
          if (is.finite(max_te_ma) && max_te_ma > 0) {
            te_ma <- 100 * te_ma / (len/bin_gene)
            plot(te_hist$mids, te_ma, type="b", col="#0066aa", lwd=4, ylim=c(0, max(te_ma)), yaxt="n", xlab="", ylab="")
            # plot(te_hist$mids, te_ma, type="b", col="#0066aa", lwd=4, ylim=c(0, max_te_ma), yaxt="n", xlab="", ylab="")
            axis(4, col="#0066aa", line = 0, col.axis = "#0066aa")
            mtext("TE+REP dens in 50 bins", side = 2, line = 2, col = "#0066aa", cex = cex_factor* 0.5, at = max(te_ma) * 0.5, adj = 0.5)
          }
        }
      }
    }
    
    # # HiC
    # if (lookup_and_plot_hic) {
    #   hic_f <- hic_files[grepl(gsub("\\|","_", chr), hic_files)]
    #   if (length(hic_f)) {
    #     hic <- read.table(hic_f, header = TRUE, sep = "\t")
    #     par(new = TRUE)
    #     plot(hic$Bin_Midpoint_BP, hic$Std_Dev_Interchrom_Contacts,
    #          type="b", col="#bb3300", lwd=4, yaxt="n", xlab="", ylab="")
    #     mtext("      HiC normalised signal", side = 2, line = 4, col = "#bb3300", cex = cex_factor* 0.5, at = 50, adj = 0)
    #   }
    # }
    
    # Gene valley
    if (!no_heli && nrow(gen_chr)) {
      gen_coords <- unlist(mapply(`:`, gen_chr$start, gen_chr$end))
      gen_coords <- gen_coords[gen_coords <= len & gen_coords >= 1]
      gen_coords <- unique(gen_coords)
      if (length(gen_coords) > 0) {
        gen_hist <- hist(gen_coords, breaks = seq(0, len, length.out = bin_gene), plot = FALSE)
        # Only apply moving average if there are enough data points (need at least 5 for padding)
        if (length(gen_hist$counts) >= 5) {
          gen_ma   <- ma(c(gen_hist$counts[1], gen_hist$counts[1], gen_hist$counts,
                           gen_hist$counts[length(gen_hist$counts)], gen_hist$counts[length(gen_hist$counts)]))[3:(length(gen_hist$counts)+2)]
          # Save gene density
          chr_density_data$gene_density <- gen_ma
          chr_density_data$gene_density_mids <- gen_hist$mids
          par(new = TRUE)
          max_gen_ma <- max(gen_ma, na.rm = TRUE)
          if (is.finite(max_gen_ma) && max_gen_ma > 0) {
            gen_ma <- 100 * gen_ma / (len/bin_gene)
            plot(gen_hist$mids, gen_ma, type="b", col="#00bb33", lwd=4, ylim=c(0, max(gen_ma)), yaxt="n", xlab="", ylab="")
            # plot(gen_hist$mids, gen_ma, type="b", col="#00bb33", lwd=4, ylim=c(0, max_gen_ma), yaxt="n", xlab="", ylab="")
            axis(4, col="#00bb33", line = 2, col.axis = "#00bb33")
            mtext("GENE dens in 50 bins", side = 2, line = 3, col = "#00bb33", cex = cex_factor* 0.5, at = max(gen_ma) * 0.5, adj = 0.5)
          }
        }
      }
    }
    
    # Store chromosome density data
    chromosome_densities[[chr]] <- chr_density_data
  }
  
  dev.off()
  print(paste0("DONE: ", plot_name))
  
  
  
  
  
}

write.csv(classes_to_plot, file = file.path(paste0(assembly_name, "_CAP_repeat_families.csv")), row.names = FALSE)

# ------------------------------------------------------------------ #
#  FINALIZE PLOT DATA WITH COMPUTED DENSITIES
# ------------------------------------------------------------------ #

plot_data <- list(
  assembly_name = assembly_name,
  chromosomes_sets = chromosomes_sets,
  chromosomes_len_sets = chromosomes_len_sets,
  scores = scores,
  gc_data = gc_data,
  ctw_data = ctw_data,
  classes_to_plot = classes_to_plot,
  edta_classes = edta_classes,
  edta_classes_colours = edta_classes_colours,
  chromosome_densities = chromosome_densities,
  no_edta = no_edta,
  no_heli = no_heli
)


# ------------------------------------------------------------------ #
#  Make the text output data
# ------------------------------------------------------------------ #
chromosomes     <- chr_info$chromosome.name
chromosomes_len <- chr_info$size
chromosome_CAP_data <- data.frame(chromosome = chromosomes, 
                                  chromosome_length = chromosomes_len)
chromosome_CAP_data$top_scoring_class = ""
chromosome_CAP_data$probability_of_top_scoring_class_is_centromeric = 0
for(j in seq_along(chromosomes)) {
  chr_classes <- scores[scores$chromosome == chromosomes[j], ]
  if(nrow(chr_classes) == 0) next
  chromosome_CAP_data$top_scoring_class[j] <- chr_classes$class[which.max(chr_classes$probability_centromeric)]
  chromosome_CAP_data$probability_of_top_scoring_class_is_centromeric[j] <- round(max(chr_classes$probability_centromeric))
}
cen_sats <- unique(scores$class[scores$probability_centromeric > 50])

summary_line1 <- paste0(assembly_name, ": probability of monocentric satellite organisation: ", round(mean(chromosome_CAP_data$probability_of_top_scoring_class_is_centromeric)), "%")
summary_line2 <- paste0("Chromosomes with over 50% probability of monocentric satellite: ", 
                        sum(chromosome_CAP_data$probability_of_top_scoring_class_is_centromeric > 50), "/", 
                        nrow(chromosome_CAP_data))
summary_line3 <- paste0("Most likely centromeric families: ")

# Prepare all output lines
output_file <- file.path(paste0(assembly_name, "_CAP_model.txt"))
output_lines <- c(summary_line1, summary_line2, summary_line3)

for(j in seq_along(cen_sats)) {
  summary_line_extra <- paste0(cen_sats[j], ": present on ", length(unique(scores$chromosome[scores$class == cen_sats[j]])), 
                               " chromosomes, with average score of ", round(mean(scores$probability_centromeric[scores$class == cen_sats[j]])),
                               "%, incliding ", sum(scores$probability_centromeric[scores$class == cen_sats[j]] > 50), " chromosomes with score above 50%")
  output_lines <- c(output_lines, summary_line_extra)
}
output_lines <- c(output_lines, "########################", "")

# Write all lines at once
writeLines(output_lines, output_file)

# Append the data frame
suppressWarnings(write.table(chromosome_CAP_data, file = output_file, append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE))

# save the data
saveRDS(plot_data, file = file.path(paste0(assembly_name, "_CAP_Rdata.rds")))
# ------------------------------------------------------------------ #

