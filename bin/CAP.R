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
    edta$start <- edta$start + 1; edta$end <- edta$end + 1
  }
  edta$start[edta$start == 0] <- 1
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
message("Removed from plotting chromosomes shorter than 200 kbp: ", chr_info$chromosome.name[chr_info$size < 200000])
chr_info <- chr_info[chr_info$size > 200000,] # TODO: remove this, user should provide fasta with only relevant chromosomes or add parameter
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
cat("Classes to plot:", classes_to_plot, "\n")





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

for(k in 1 : length(chromosomes_sets)) {
  
  plot_name <- file.path(paste0(assembly_name, "_CAP_plot_", k, "_", suffix, ".png"))
  
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
  cex_factor <- 2.2
  # ------------------------------------------------------------------ #
  #  TITLE
  # ------------------------------------------------------------------ #
  plot(NA, xlim = c(1,100), ylim = c(1,100), axes = FALSE, xlab = "", ylab = "")
  text(50, 95, assembly_name, pos = 1, cex = cex_factor*5)
  text(1, 70, sprintf("Total size: %s Mbp", formatC(sum(chromosomes_len)/1e6, format="f", digits=3, big.mark=" ")), pos=4, cex = cex_factor*3)
  text(1,55, sprintf("Chromosomes: %d", length(chromosomes)), pos=4, cex = cex_factor*3)
  text(1,40, sprintf("Transposable elements legend:"), pos=4, cex = cex_factor*2)
  
  text(x = c(1,15,30,45,60,75,90), y = 32, 
       labels =  c("class I LTR: ",
                   "Gypsy",
                   "Copia",
                   "Bel Pao",
                   "TRIM",
                   "Caulimoviridae", 
                   "unspecified"), 
       col = c("black", edta_classes_colours[1:6]),
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15,30,45,60,75), y = 26, 
       labels = c("class I non-LTR: ", 
                  "LINE",
                  "SINE",
                  "Penelope",
                  "DIRS YR",
                  "unspecified"),
       col = c("black", edta_classes_colours[7:11]),
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15), y = 20, 
       labels = c("class II TIRs: ", "Kolobok;Ginger;Academ;Novosib;Sola;Merlin;IS3EU;PiggyBac;hAT;Mutator;Tc1 Mariner;Dada;CACTA;Zisupton;PIF Harbinger"),
       col = c("black", edta_classes_colours[12]), 
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15,30,45,60), y = 14, 
       labels = c("class II others: ",  
                  "DNA_transposon",
                  "helitron",
                  "MITE",
                  "Maverick Polinton"),
       col = c("black", edta_classes_colours[13:16]), 
       pos = 4, cex = cex_factor* 1.5)
  text(x = c(1,15,30), y = 8, 
       labels = c("other: ", "Tyrosine Recombinase",
                  "unspecified"),
       col = c("black", edta_classes_colours[17:18]), 
       pos = 4, cex = cex_factor* 1.5)
  
  # ------------------------------------------------------------------ #
  #  PER-CHROMOSOME LOOP
  # ------------------------------------------------------------------ #
  
  for (j in seq_along(chromosomes)) {
    chr <- chromosomes[j]
    len <- chromosomes_len[j]
    rep_chr <- subset(repeats, seqID == chr)
    edt_chr <- if (!no_edta) subset(edta, seqID == chr) else data.frame()
    gen_chr <- if (!no_heli) subset(genes, seqID == chr) else data.frame()
    
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
    } else plot.new()
    
    # === PLOT A: Repeats + GC + families ===
    win_rep <- genomic.bins.starts(1, len, bin.size = bin_rep)
    rep_cov <- if (nrow(rep_chr)) {
      calculate.repeats.percentage.in.windows(win_rep, rep_chr$start, rep_chr$width, len)
    } else rep(0, length(win_rep))
    rep_cov[rep_cov == 0] <- NA; rep_cov[rep_cov > 100] <- 100
    plot(win_rep, rep_cov, type="h", col="#CCCCCC", ylim=c(0,100), xaxt="n", yaxt="n", xlab="", ylab="")
    mtext("REP%         per 10 Kbp", side = 2, line = 0, col = "grey", cex = cex_factor* 0.5, at = 10, adj = 0)
    axis(side = 2, labels = c("0","100"), at = c(0,100))
    
    # GC (TODO only if requested)
    gc_chs_data <- gc_data[gc_data$chromosome == chr,]
    gc_mids <- gc_chs_data$bin_mid
    gc_vals <- gc_chs_data$bin_value
    lines(gc_mids, gc_vals, type = "l", lwd = 2, cex = cex_factor * 1)
    mtext("GC%            per  2 Kbp", side = 2, line = 2, col = "black", cex = cex_factor* 0.5, at = 10, adj = 0)
    
    # CTW (TODO only if requested)
    ctw_chs_data <- ctw_data[ctw_data$chromosome == chr,]
    ctw_mids <- ctw_chs_data$bin_mid
    ctw_vals <- ctw_chs_data$bin_value
    lines(x = ctw_mids, ctw_vals, type = "l", col = "#FFA500", lwd = 2, cex = cex_factor * 1)
    mtext("CTW            per  2 Kbp D10", side = 2, line = 3, col = "#FFA500", cex = cex_factor* 0.5, at = 10, adj = 0)
    
    
    # Families
    if (length(classes_to_plot)) {
      for (k in seq_len(length(classes_to_plot))) {
        fam <- subset(rep_chr, new_class == classes_to_plot[k])
        if (!nrow(fam)) next
        cov <- calculate.repeats.percentage.in.windows(win_rep, fam$start, fam$width, len)
        cov[cov == 0] <- NA; cov[cov > 100] <- 100
        lines(win_rep, cov, col = palette[k], pch=16, type="o", lwd = 2, cex = cex_factor * 2)
      }
      mtext("SIG REP% per 10 Kbp", side = 2, line = 1, col = "#88CCEE", cex = cex_factor* 0.5, at = 10, adj = 0)
    }
    axis(1, at = seq(0, len, by = x_tick), labels = FALSE)
    # text(1, bin_rep, chr, pos = 4)
    
    # === PLOT B: EDTA + TE/repeat peak + genes + HiC ===
    win_edta <- genomic.bins.starts(1, len, bin.size = bin_edta)
    edt_cov <- if (!no_edta && nrow(edt_chr)) {
      calculate.repeats.percentage.in.windows(win_edta, edt_chr$start, edt_chr$width, len)
    } else rep(0, length(win_edta))
    edt_cov[edt_cov == 0] <- NA; edt_cov[edt_cov > 100] <- 100
    plot(win_edta + bin_edta/2, edt_cov, type="h", col="#CCCCCC", ylim=c(0,100), xaxt="n", yaxt="n", xlab="", ylab="")
    mtext("EDTA%           per 100 Kbp", side = 2, line = 0, col = "grey", cex = cex_factor* 0.5, at = 10, adj = 0)
    axis(2, col="grey", at = c(0,100), labels = c("0", "100"))
    
    # EDTA classes
    if (!no_edta && nrow(edt_chr)) {
      for (k in rev(seq_along(edta_classes))) {
        cls <- edt_chr[edt_chr$type %in% edta_classes[[k]], ]
        if (!nrow(cls)) next
        cov <- calculate.repeats.percentage.in.windows(win_edta, cls$end, cls$width, len)
        cov[cov == 0] <- NA; cov[cov > 100] <- 100
        lines(win_edta + (bin_edta/2), cov, col = edta_classes_colours[k], pch=16, type="o", lwd = 2, cex = cex_factor * 1)
      }
      mtext("FAM EDTA%  per 100 Kbp", side = 2, line = 1, col = "red", cex = cex_factor* 0.5, at = 10, adj = 0)
    }
    
    # TE+repeat peak
    te_coords <- c()
    if (!no_edta && nrow(edt_chr)) te_coords <- c(te_coords, unlist(mapply(`:`, edt_chr$start, edt_chr$end)))
    if (nrow(rep_chr))          te_coords <- c(te_coords, unlist(mapply(`:`, rep_chr$start, rep_chr$end)))
    if (length(te_coords)) {
      te_coords <- te_coords[te_coords <= len & te_coords >= 1]
      te_hist <- hist(te_coords, breaks = seq(0, len, length.out = bin_gene), plot = FALSE)
      te_ma   <- ma(c(te_hist$counts[1], te_hist$counts[1], te_hist$counts,
                      te_hist$counts[length(te_hist$counts)], te_hist$counts[length(te_hist$counts)]))[3:(length(te_hist$counts)+2)]
      par(new = TRUE)
      plot(te_hist$mids, te_ma, type="b", col="#0066aa", lwd=4, ylim=c(0, max(te_ma)), yaxt="n", xlab="", ylab="")
      axis(4, col="#0066aa", line = 0, col.axis = "#0066aa")
      mtext("      TE+REP dens per 100 Kbp", side = 2, line = 2, col = "#0066aa", cex = cex_factor* 0.5, at = 20, adj = 0)
    }
    
    # # HiC
    # if (lookup_and_plot_hic) {
    #   hic_f <- hic_files[grepl(gsub("\\|","_", chr), hic_files)]
    #   if (length(hic_f)) {
    #     hic <- read.table(hic_f, header = TRUE, sep = "\t")
    #     par(new = TRUE)
    #     plot(hic$Bin_Midpoint_BP, hic$Std_Dev_Interchrom_Contacts,
    #          type="b", col="#bb3300", lwd=4, yaxt="n", xlab="", ylab="")
    #     mtext("      HiC normalised signal", side = 2, line = 4, col = "#bb3300", cex = cex_factor* 0.5, at = 20, adj = 0)
    #   }
    # }
    
    # Gene valley
    if (!no_heli && nrow(gen_chr)) {
      gen_coords <- unlist(mapply(`:`, gen_chr$start, gen_chr$end))
      gen_coords <- gen_coords[gen_coords <= len & gen_coords >= 1]
      gen_hist <- hist(gen_coords, breaks = seq(0, len, length.out = bin_gene), plot = FALSE)
      gen_ma   <- ma(c(gen_hist$counts[1], gen_hist$counts[1], gen_hist$counts,
                       gen_hist$counts[length(gen_hist$counts)], gen_hist$counts[length(gen_hist$counts)]))[3:(length(gen_hist$counts)+2)]
      par(new = TRUE)
      plot(gen_hist$mids, gen_ma, type="b", col="#00bb33", lwd=4, ylim=c(0, max(gen_ma)), yaxt="n", xlab="", ylab="")
      axis(4, col="#00bb33", line = 2, col.axis = "#00bb33")
      mtext("      GENE dens     per 100 Kbp", side = 2, line = 3, col = "#00bb33", cex = cex_factor* 0.5, at = 10, adj = 0)
    }
  }
  
  dev.off()
  message("DONE: ", plot_name)
  
  
  
  
  
}


# ------------------------------------------------------------------ #
#  DOT-PLOT & FINAL TABLE
# ------------------------------------------------------------------ #

if(T) {
  classes_to_plot <- data.frame(class = classes_to_plot,
                                consensus = classes$consensus[match(classes_to_plot, classes$class)],
                                count = classes$count[match(classes_to_plot, classes$class)],
                                stringsAsFactors = FALSE)
  
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
            cex = cex_factor* 10, xaxt = "n", yaxt = "n", lwd = 3)
    
    abline(v = c(1, divs, divs_rc), col = classes_to_plot$colour, lwd = 12)
    abline(h = c(1, divs, divs_rc), col = classes_to_plot$colour, lwd = 12)
    abline(v = rev(divs_rc)[1:nrow(classes_to_plot)], col = classes_to_plot$colour, lwd = 12)
    abline(h = rev(divs_rc)[1:nrow(classes_to_plot)], col = classes_to_plot$colour, lwd = 12)
    abline(v = nchar(cons_seq), h = nchar(cons_seq), col = "black", lwd = 24)
    
    dev.off()
    message("DONE: ", plot_name)
    
  }
  
}


write.csv(classes_to_plot, file = file.path(paste0(assembly_name, "_CAP_repeat_families.csv")), row.names = FALSE)
write.csv(classes_to_plot, file = file.path(paste0(assembly_name, "_CAP_model.txt")), row.names = FALSE)



















# ------------------------------------------------------------------ #
#  Make the text output data
# ------------------------------------------------------------------ #





















# Save outputs
# write.csv(repeat_families_data, file = output_cap_repeat_families, row.names = FALSE)  # Assuming created in code
# writeLines(model_text, con = output_cap_model)  # Assuming text output