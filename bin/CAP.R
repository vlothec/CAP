#!/usr/bin/env Rscript

# CAP.R
# Description: Generate final CAP outputs with optional TEs and genes

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9 || length(args) > 11) {
  stop("Usage: CAP.R <predictions.csv> <repeats_reclassed.csv> <arrays_reclassed.csv> <genome_classes.csv> <metadata.csv> <assembly_name> <GC> <CTW> [<TEs_filtered.csv>] [<genes_filtered.csv>]")
}
no_edta <- FALSE
no_heli <- FALSE

predictions_csv <- args[1]
repeats_reclassed_csv <- args[2]
arrays_reclassed_csv <- args[3]
genome_classes_csv <- args[4]
metadata_csv <- args[5]
assembly_name <- basename(args[6])
gc_csv <- basename(args[7])
ctw_csv <- basename(args[8])
tes_filtered_csv <- if (args[9] != "NO_FILE") args[9] else no_edta <- TRUE
genes_filtered_csv <- if (args[10] != "NO_FILE") args[10] else no_heli <- TRUE
scores_csv <- basename(args[11])


# Load libraries
suppressMessages({
  library(seqinr)
  library(stringr)
  library(Biostrings)
  library(GenomicRanges)
  library(scales)
  library(BCT)
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

# Load optionals if provided
print("load extra")
if (!no_edta) tes_data <- read.csv(tes_filtered_csv)
if (!no_heli) genes_data <- read.csv(genes_filtered_csv)







repeats  <- repeats_data
arrays   <- arrays_data
classes  <- genome_classes_data
chr_info <- metadata_data
if (nrow(classes)) classes$num_ID <- seq_len(nrow(classes))

if (!no_heli) {
  genes <- genes_data
  genes <- genes[genes$V3 == "CDS", ]
}

if (!no_edta) {
  edta <- tes_data
  if ("reassigned" %in% colnames(edta)) {
    names(edta) <- c("", "V1","V2","V3","V4","V5","V6","V7","V8","ID","Name",
                     "Classification","Sequence_ontology","Identity","Method",
                     "TSD","TIR","motif","tsd","oldV3","overlapping_bp",
                     "width","overlapping_percentage","reassigned")
    edta$V4 <- edta$V4 + 1; edta$V5 <- edta$V5 + 1
  }
  edta$V4[edta$V4 == 0] <- 1
}

####
edta_classes <- list(
  # class I (retrotransposons)
  ## LTR retrotransposons
  c("Gypsy_LTR_retrotransposon"),
  c("Copia_LTR_retrotransposon"),
  c("Bel_Pao_LTR_retrotransposon"),
  c("TRIM_LTR_retrotransposon"),
  c("Caulimoviridae"),
  c("Retrovirus", "LTR_retrotransposon", "long_terminal_repeat"),
  ## Non-LTR retrotransposons
  c("LINE_element"),
  c("SINE_element"),
  c("Penelope_retrotransposon"),
  c("DIRS_YR_retrotransposon"),
  c("non_LTR_retrotransposon"),
  # class II (DNA transposons)
  ## TIRs
  c("Kolobok_TIR_transposon" , "Ginger_TIR_transposon", "Academ_TIR_transposon", "Novosib_TIR_transposon", "Sola_TIR_transposon", "Merlin_TIR_transposon", "IS3EU_TIR_transposon", "PiggyBac_TIR_transposon", "hAT_TIR_transposon", "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon", "Dada_TIR_transposon", "CACTA_TIR_transposon", "Zisupton_TIR_transposon", "PIF_Harbinger_TIR_transposon"),
  ## other class II
  c("DNA_transposon"),
  c("helitron"),
  c("MITE"),
  c("Maverick_Polinton", "polinton"),
  # other, recombinase element based
  c("Tyrosine_Recombinase_Elements", "Crypton_Tyrosine_Recombinase"),
  # others
  c("TE", "TE_unclass"),
  # likely not TEs, remove for plotting?
  c("repeat_region", "SUPER", "Sequence_Ontology", "rRNA_gene", "target_site_duplication", "chr"))
edta_classes_colours <-  c(
  "#E31A1C",  # red
  "#D55E00",  # reddish-orange
  "#F5793A",  # bright orange
  "#FF7F00",  # orange
  "#FDBF6F",  # peach
  "#F0E442",  # yellow
  "#6A3D9A",  # dark purple
  "#A95AA1",  # purple
  "#CC79A7",  # pink
  "#DDA0DD",  # light pinkish purple (plum)
  "#CAB2D6",  # lavender
  "#1F78B4",  # blue
  "#B2DF8A",  # light green
  "#33A02C",  # green
  "#009E73",  # teal green
  "#56B4E9",  # light blue
  "#7F7F7F",  # grey
  "#000000",   # black
  "#000000"   # black
)

# ------------------------------------------------------------------ #
#  CHROMOSOME FILTERING (top 30 longest)
# ------------------------------------------------------------------ #
chromosomes     <- chr_info$chromosome.name
chromosomes_len <- chr_info$size

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
scores <- cbind(scores_data, predictions_data)

if (nrow(scores)) {
  scores <- scores[order(scores$chromosome), ]
  # ---- rescoring (unchanged logic) ---------------------------------
  rescoring <- function(df) {
    size_norm <- df$total_bp / max(df$total_bp)
    size_norm <- size_norm + (1 - size_norm) / 1.5
    df$score_total <- 0
    
    add <- function(col, cap, w) {
      v <- pmin(df[[col]], cap); v[v < 0] <- 0
      df$score_total <<- df$score_total + (v / cap) * w * size_norm
    }
    add("total_bp_norm_chr",        0.1, 2)
    add("total_bp_norm_rep",        0.5, 1)
    add("start_sd_norm_chr",        0.15, 1)
    add("start_norm_chr_0_50",      0.25, 1)
    add("gaps_with_TEs_fraction",   0.75, 0.5)
    add("centre_array_edit",        20, 2)
    add("centre_array_width_sd",    20, 2)
    add("centre_chromosome_edit",   15, 1)
    add("centre_chromosome_width_sd",15, 1)
    
    df$TE_prox_score   <- with(df, abs(TE_lm_coef) / (1 + TE_prox_dist + TE_prox_SD) * 100)
    df$gene_prox_score <- with(df, abs(gene_lm_coef) / (1 + gene_prox_dist + gene_prox_SD) * 100)
    add("TE_prox_score",   5, 5)
    add("gene_prox_score", 5, 5)
    
    df
  }
  scores <- do.call(rbind, by(scores, scores$chromosome, rescoring))
  scores <- scores[, -1]
  
  scores$ed_perc <- 100 * scores$centre_array_edit / scores$mean_length
  scores$width_sd_perc <- 100 * scores$centre_array_width_sd / scores$mean_length
  scores$ed_perc[scores$ed_perc < 0] = NA
  scores$width_sd_perc[scores$width_sd_perc < 0] = NA
  scores$ed_perc = 100 - scores$ed_perc
  scores$width_sd_perc = 100 - scores$width_sd_perc
}

# ------------------------------------------------------------------ #
#  SELECT CLASSES TO PLOT (top-scoring, >=5 kbp, <=20)
# ------------------------------------------------------------------ #

classes_to_plot <- if (nrow(scores)) {
  filter_top_classes <- function(s, c) {
    top <- do.call(rbind, by(s, s$chromosome, function(df) {
      df <- df[order(df$score_total, decreasing = TRUE), ]
      df <- df[df$score_total >= 3, ]
      if (nrow(df) > 4) df <- df[1:4, ]
      df$total_genome_bp <- c$sum_coverage[match(df$class, c$class)]
      df[df$total_genome_bp >= 5000, ]
    }))
    c$to_plot <- c$class %in% top$class
    c <- c[c$to_plot, ]
    if (nrow(c) > 20) c <- c[1:20, ]
    c
  }
  filter_top_classes(scores, classes)
} else classes[FALSE, ]



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

for(k in 1 : length(chromosomes_sets)) {
  
  plot_name <- file.path(paste0(assembly_name, "_CAP_plot_", k, "_", suffix, ".png"))
  
  
  # if (file.exists(plot_name)) {
  #   if(rerun_if_completed) {
  #     message("Plot exists: ", plot_name, ", removing and plotting again")
  #     file.remove(plot_name)
  #   } else {
  #     message("Plot exists: ", plot_name, ", moving on, use --rerun to remove existing and plot again")
  #     next
  #   }
  #   
  # }
  
  
  chromosomes <- chromosomes_sets[[k]]
  chromosomes_len <- chromosomes_len_sets[[k]]
  
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
  
  # ------------------------------------------------------------------ #
  #  TITLE
  # ------------------------------------------------------------------ #
  plot(NA, xlim = c(1,100), ylim = c(1,100), axes = FALSE, xlab = "", ylab = "")
  text(50, 95, assembly_name, pos = 1, cex = 5)
  text(1, 70, sprintf("Total size: %s Mbp", formatC(sum(chromosomes_len)/1e6, format="f", digits=3, big.mark=" ")), pos=4, cex=3)
  text(1,55, sprintf("Chromosomes: %d", length(chromosomes)), pos=4, cex=3)
  text(1,40, sprintf("Transposable elements legend:"), pos=4, cex=2)
  
  text(x = c(1,15,30,45,60,75,90), y = 32, 
       labels =  c("class I LTR: ",
                   "Gypsy",
                   "Copia",
                   "Bel Pao",
                   "TRIM",
                   "Caulimoviridae", 
                   "unspecified"), 
       col = c("black", edta_classes_colours[1:6]),
       pos = 4, cex = 1.5)
  text(x = c(1,15,30,45,60,75), y = 26, 
       labels = c("class I non-LTR: ", 
                  "LINE",
                  "SINE",
                  "Penelope",
                  "DIRS YR",
                  "unspecified"),
       col = c("black", edta_classes_colours[7:11]),
       pos = 4, cex = 1.5)
  text(x = c(1,15), y = 20, 
       labels = c("class II TIRs: ", "Kolobok;Ginger;Academ;Novosib;Sola;Merlin;IS3EU;PiggyBac;hAT;Mutator;Tc1 Mariner;Dada;CACTA;Zisupton;PIF Harbinger"),
       col = c("black", edta_classes_colours[12]), 
       pos = 4, cex = 1.5)
  text(x = c(1,15,30,45,60), y = 14, 
       labels = c("class II others: ",  
                  "DNA_transposon",
                  "helitron",
                  "MITE",
                  "Maverick Polinton"),
       col = c("black", edta_classes_colours[13:16]), 
       pos = 4, cex = 1.5)
  text(x = c(1,15,30), y = 8, 
       labels = c("other: ", "Tyrosine Recombinase",
                  "unspecified"),
       col = c("black", edta_classes_colours[17:18]), 
       pos = 4, cex = 1.5)
  
  # ------------------------------------------------------------------ #
  #  PER-CHROMOSOME LOOP
  # ------------------------------------------------------------------ #
  for (j in seq_along(chromosomes)) {
    chr <- chromosomes[j]; len <- chromosomes_len[j]
    rep_chr <- subset(repeats, seqID == chr)
    edt_chr <- if (!no_edta) subset(edta, V1 == chr) else data.frame()
    gen_chr <- if (!no_heli) subset(genes, V1 == chr) else data.frame()
    
    print(paste0("Genome ", assembly_name, " | Chromosome ", j, "/", length(chromosomes)))
    
    plot(NA,NA, xlim = c(1,100), ylim = c(1,300), xlab = "", ylab = "", axes = F)
    text(x = 1, y = 15, labels = assembly_name, pos = 4, cex = 1.4)
    text(x = 15, y = 15, labels = chromosomes[j], pos = 4, cex = 1.2)
    text(x = 30, y = 15, 
         labels = paste0(formatC(chromosomes_len[j]/1000000, format = "f", big.mark = " ", digits = 3), " Mbp"), 
         pos = 4, cex = 1.2)
    abline(h = 50, lwd = 4)
    abline(h = 60, lwd = 4)
    
    
    # === TABLE ===
    if (nrow(scores)) {
      sc <- subset(scores, chromosome == chr)[, c("class","count","mean_length","total_bp",
                                                  "ed_perc","width_sd_perc","score_total")]
      sc[, 3:7] <- round(sc[, 3:7], 2)
      sc$colours <- palette[match(sc$class, classes_to_plot$class)]
      create_table(sc[,1:7], c("Class","Repeats no","Mean width, bp","Total bp",
                               "Sequence similarity %","Width similarity %","Centromeric score"),
                   colours = sc$colours)
    } else plot.new()
    
    # === PLOT A: Repeats + GC + families ===
    win_rep <- genomic.bins.starts(1, len, bin.size = bin_rep)
    rep_cov <- if (nrow(rep_chr)) {
      calculate.repeats.percentage.in.windows(win_rep, rep_chr$start, rep_chr$width, len)
    } else rep(0, length(win_rep))
    rep_cov[rep_cov == 0] <- NA; rep_cov[rep_cov > 100] <- 100
    plot(win_rep, rep_cov, type="h", col="#CCCCCC", ylim=c(0,100), xaxt="n", yaxt="n", xlab="", ylab="")
    mtext("REP%         per 10 Kbp", side = 2, line = 0, col = "grey", cex = 0.5, at = 10, adj = 0)
    axis(side = 2, labels = c("0","100"), at = c(0,100))
    
    # GC (TODO only if requested)
    gc_chs_data <- gc_data[gc_data$chromosome == chr,]
    gc_mids <- gc_chs_data$bin_mid
    gc_vals <- gc_chs_data$bin_value
    lines(gc_mids, gc_vals, type = "l", lwd = 1)
    mtext("GC%            per  2 Kbp", side = 2, line = 2, col = "black", cex = 0.5, at = 10, adj = 0)
    
    # CTW (TODO only if requested)
    ctw_chs_data <- ctw_data[ctw_data$chromosome == chr,]
    ctw_mids <- ctw_chs_data$bin_mid
    ctw_vals <- ctw_chs_data$bin_value
    lines(x = ctw_mids, ctw_vals, type = "l", col = "#FFA500",lwd = 1)
    mtext("CTW            per  2 Kbp D10", side = 2, line = 3, col = "#FFA500", cex = 0.5, at = 10, adj = 0)
    
    
    # Families
    if (nrow(classes_to_plot)) {
      for (k in seq_len(nrow(classes_to_plot))) {
        fam <- subset(rep_chr, new_class == classes_to_plot$class[k])
        if (!nrow(fam)) next
        cov <- calculate.repeats.percentage.in.windows(win_rep, fam$start, fam$width, len)
        cov[cov == 0] <- NA; cov[cov > 100] <- 100
        lines(win_rep, cov, col = palette[k], pch=16, type="o")
      }
      mtext("SIG REP% per 10 Kbp", side = 2, line = 1, col = "#88CCEE", cex = 0.5, at = 10, adj = 0)
    }
    axis(1, at = seq(0, len, by = x_tick), labels = FALSE)
    # text(1, bin_rep, chr, pos = 4)
    
    # === PLOT B: EDTA + TE/repeat peak + genes + HiC ===
    win_edta <- genomic.bins.starts(1, len, bin.size = bin_edta)
    edt_cov <- if (!no_edta && nrow(edt_chr)) {
      calculate.repeats.percentage.in.windows(win_edta, edt_chr$V4, edt_chr$width, len)
    } else rep(0, length(win_edta))
    edt_cov[edt_cov == 0] <- NA; edt_cov[edt_cov > 100] <- 100
    plot(win_edta + bin_edta/2, edt_cov, type="h", col="#CCCCCC", ylim=c(0,100), xaxt="n", yaxt="n", xlab="", ylab="")
    mtext("EDTA%           per 100 Kbp", side = 2, line = 0, col = "grey", cex = 0.5, at = 10, adj = 0)
    axis(2, col="grey", at = c(0,100), labels = c("0", "100"))
    
    # EDTA classes
    if (!no_edta && nrow(edt_chr)) {
      for (k in rev(seq_along(edta_classes))) {
        cls <- edt_chr[edt_chr$V3 %in% edta_classes[[k]], ]
        if (!nrow(cls)) next
        cov <- calculate.repeats.percentage.in.windows(win_edta, cls$V4, cls$width, len)
        cov[cov == 0] <- NA; cov[cov > 100] <- 100
        lines(win_edta + (bin_edta/2), cov, col = edta_classes_colours[k], pch=16, type="o")
      }
      mtext("FAM EDTA%  per 100 Kbp", side = 2, line = 1, col = "red", cex = 0.5, at = 10, adj = 0)
    }
    
    # TE+repeat peak
    te_coords <- c()
    if (!no_edta && nrow(edt_chr)) te_coords <- c(te_coords, unlist(mapply(`:`, edt_chr$V4, edt_chr$V5)))
    if (nrow(rep_chr))          te_coords <- c(te_coords, unlist(mapply(`:`, rep_chr$start, rep_chr$end)))
    if (length(te_coords)) {
      te_hist <- hist(te_coords, breaks = seq(0, len, length.out = bin_gene), plot = FALSE)
      te_ma   <- ma(c(te_hist$counts[1], te_hist$counts[1], te_hist$counts,
                      te_hist$counts[length(te_hist$counts)], te_hist$counts[length(te_hist$counts)]))[3:(length(te_hist$counts)+2)]
      par(new = TRUE)
      plot(te_hist$mids, te_ma, type="b", col="#0066aa", lwd=4, ylim=c(0, max(te_ma)), yaxt="n", xlab="", ylab="")
      axis(4, col="#0066aa", line = 0, col.axis = "#0066aa")
      mtext("      TE+REP dens per 100 Kbp", side = 2, line = 2, col = "#0066aa", cex = 0.5, at = 20, adj = 0)
    }
    
    # HiC
    if (lookup_and_plot_hic) {
      hic_f <- hic_files[grepl(gsub("\\|","_", chr), hic_files)]
      if (length(hic_f)) {
        hic <- read.table(hic_f, header = TRUE, sep = "\t")
        par(new = TRUE)
        plot(hic$Bin_Midpoint_BP, hic$Std_Dev_Interchrom_Contacts,
             type="b", col="#bb3300", lwd=4, yaxt="n", xlab="", ylab="")
        mtext("      HiC normalised signal", side = 2, line = 4, col = "#bb3300", cex = 0.5, at = 20, adj = 0)
      }
    }
    
    # Gene valley
    if (!no_heli && nrow(gen_chr)) {
      gen_coords <- unlist(mapply(`:`, gen_chr$V4, gen_chr$V5))
      gen_hist <- hist(gen_coords, breaks = seq(0, len, length.out = bin_gene), plot = FALSE)
      gen_ma   <- ma(c(gen_hist$counts[1], gen_hist$counts[1], gen_hist$counts,
                       gen_hist$counts[length(gen_hist$counts)], gen_hist$counts[length(gen_hist$counts)]))[3:(length(gen_hist$counts)+2)]
      par(new = TRUE)
      plot(gen_hist$mids, gen_ma, type="b", col="#00bb33", lwd=4, ylim=c(0, max(gen_ma)), yaxt="n", xlab="", ylab="")
      axis(4, col="#00bb33", line = 2, col.axis = "#00bb33")
      mtext("      GENE dens     per 100 Kbp", side = 2, line = 3, col = "#00bb33", cex = 0.5, at = 10, adj = 0)
    }
  }
  
  dev.off()
  message("DONE: ", plot_name)
  
  
  
  
  
}


# ------------------------------------------------------------------ #
#  DOT-PLOT & FINAL TABLE
# ------------------------------------------------------------------ #
if(T) {
  
  # if (nrow(classes_to_plot)) {
  #   # plot_name <- output_cap_plot
  #   plot_name <- file.path(paste0("repeat_total_table_", assembly_name, ".png"))
  #   png(plot_name, width = 3000, height = 200 + 100 * nrow(classes_to_plot))
  #   create_table(classes_to_plot[, c("class","count","median_length","sum_coverage")],
  #                c("Class","Repeats no","Median width", "Total width"),
  #                colours = palette[1:nrow(classes_to_plot)], font_size = 4)
  #   dev.off()
  #   message("DONE: ", plot_name)
  # }
  classes_to_plot$colour <- palette[1:nrow(classes_to_plot)]
  classes_to_plot <- classes_to_plot[classes_to_plot$count > 50, ]
  
  if (nrow(classes_to_plot)) {
    
    cons_seq <- paste(classes_to_plot$consensus, collapse = "")
    rev_seq  <- revCompString(cons_seq)
    full_seq <- paste(cons_seq, rev_seq, sep = "")
    divs     <- cumsum(nchar(classes_to_plot$consensus))
    divs_rc  <- divs + nchar(cons_seq)
    
    # plot_name <- output_cap_dotplot
    plot_name <- file.path(paste0(assembly_name, "_CAP_dotplot.png"))
    png(plot_name, width = nchar(full_seq), height = nchar(full_seq))
    dotPlot(strsplit(full_seq, "")[[1]], strsplit(full_seq, "")[[1]],
            wsize = 4, wstep = 1, nmatch = 4,
            col = c("white","black"), xlab = "n", ylab = "n",
            cex = 10, xaxt = "n", yaxt = "n", lwd = 3)
    
    abline(v = c(1, divs, divs_rc), col = classes_to_plot$colour, lwd = 12)
    abline(h = c(1, divs, divs_rc), col = classes_to_plot$colour, lwd = 12)
    abline(v = rev(divs_rc)[1:nrow(classes_to_plot)], col = classes_to_plot$colour, lwd = 12)
    abline(h = rev(divs_rc)[1:nrow(classes_to_plot)], col = classes_to_plot$colour, lwd = 12)
    abline(v = nchar(cons_seq), h = nchar(cons_seq), col = "black", lwd = 24)
    
    dev.off()
    message("DONE: ", plot_name)
    
  }
  
}






















# ------------------------------------------------------------------ #
#  Make the text output data
# ------------------------------------------------------------------ #





















# Save outputs
# write.csv(repeat_families_data, file = output_cap_repeat_families, row.names = FALSE)  # Assuming created in code
# writeLines(model_text, con = output_cap_model)  # Assuming text output