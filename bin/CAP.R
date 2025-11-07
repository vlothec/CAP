#!/usr/bin/env Rscript

# CAP.R
# Description: Generate final CAP outputs with optional TEs and genes

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8 || length(args) > 10) {
  stop("Usage: CAP.R <predictions.csv> <repeats_reclassed.csv> <arrays_reclassed.csv> <genome_classes.csv> <metadata.csv> [<TEs_filtered.csv>] [<genes_filtered.csv>] <output_CAP_plot.png> <output_CAP_repeat_families.csv> <output_CAP_dotplot.png> <output_CAP_model.txt>")
}

predictions_csv <- args[1]
repeats_reclassed_csv <- args[2]
arrays_reclassed_csv <- args[3]
genome_classes_csv <- args[4]
metadata_csv <- args[5]
arg_index <- 6

tes_filtered_csv <- if (length(args) >= arg_index + 1 && args[arg_index] != "") args[arg_index] else NULL
if (!is.null(tes_filtered_csv)) arg_index <- arg_index + 1

genes_filtered_csv <- if (length(args) >= arg_index + 1 && args[arg_index] != "") args[arg_index] else NULL
if (!is.null(genes_filtered_csv)) arg_index <- arg_index + 1

output_cap_plot <- args[arg_index]
output_cap_repeat_families <- args[arg_index + 1]
output_cap_dotplot <- args[arg_index + 2]
output_cap_model <- args[arg_index + 3]

# Load libraries
suppressMessages({library(seqinr)
                  library(msa)})

# Load additional functions
source("./aux_fun.R")

# Load required data
predictions_data <- read.csv(predictions_csv)
repeats_data <- read.csv(repeats_reclassed_csv)
arrays_data <- read.csv(arrays_reclassed_csv)
genome_classes_data <- read.csv(genome_classes_csv)
metadata_data <- read.csv(metadata_csv)

# Load optionals if provided
if (!is.null(tes_filtered_csv)) {
  tes_data <- read.csv(tes_filtered_csv)
}
if (!is.null(genes_filtered_csv)) {
  genes_data <- read.csv(genes_filtered_csv)
}

# --- Your CAP generation code here (including plots) ---

# Save outputs
# For plots, use png() or ggsave()
# png(output_cap_plot); plot(...); dev.off()
write.csv(repeat_families_data, file = output_cap_repeat_families, row.names = FALSE)  # Assuming created in code
writeLines(model_text, con = output_cap_model)  # Assuming text output