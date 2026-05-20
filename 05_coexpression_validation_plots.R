# ============================================================
# Script: 05_coexpression_validation_plots.R
# Purpose: Validate selected regulator-target gene pairs by plotting
#          pairwise co-expression relationships in LIHC and GTEx liver
#          expression datasets.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files
VALIDATION_DIR <- file.path(INPUT_DIR, "validation_UCSC_LIHC_GTEx/")
RESULTS_DIR <- file.path(INPUT_DIR, "results/")

GTEX_PHENOTYPE_FILE <- file.path(VALIDATION_DIR, "GTEX_phenotype.gz")
LIHC_SURVIVAL_FILE <- file.path(VALIDATION_DIR, "survival_LIHC_survival.txt")
GENE_PROBEMAP_FILE <- file.path(VALIDATION_DIR, "probeMap_gencode.v23.annotation.gene.probemap")
LIVER_EXPRESSION_FILE <- file.path(RESULTS_DIR, "liver_tcga_gtex.rds")

# Output files
OUTPUT_COEXPRESSION_FIGURE <- file.path(OUTPUT_DIR, "selected_pair_coexpression.svg")
OUTPUT_MYC_COEXPRESSION_FIGURE <- file.path(OUTPUT_DIR, "MYC_pair_coexpression.svg")

# Plot parameters
POINT_SIZE <- 5
POINT_ALPHA <- 0.25
PLOT_DPI <- 300

# -----------------------------
# LIBRARIES
# -----------------------------
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggpmisc)

# -----------------------------
# INPUT GENE PAIRS
# -----------------------------
# Selected regulator-target pairs from the network validation analysis.
selected_gene_pairs <- list(
  c("AQR", "PES1"),
  c("RBM39", "KIF1C"),
  c("RPS10", "FASTKD1"),
  c("HNRNPC", "AKAP8"),
  c("IGF2BP1", "CCAR1"),
  c("IGF2BP1", "PCBP2"),
  c("MSI2", "PES1"),
  c("PPIL4", "LIN28B"),
  c("HSPD1", "AKAP8L"),
  c("CEBPZ", "RBM27"),
  c("EWSR1", "PES1"),
  c("AQR", "XRCC5"),
  c("EEF2", "AQR"),
  c("BOP1", "CEBPZ"),
  c("HSPD1", "PKM"),
  c("HNRNPK", "SNRNP70"),
  c("RCC2", "WRN"),
  c("HNRNPA2B1", "CIRBP"),
  c("ILF2", "CELF1")
)

# MYC-related candidate pairs.
myc_gene_pairs <- list(
  c("MYC", "ADAR"),
  c("MYC", "AKAP8L"),
  c("MYC", "CCAR1"),
  c("MYC", "EIF4G1"),
  c("MYC", "FASTKD1"),
  c("MYC", "HNRNPK"),
  c("MYC", "HNRNPLL"),
  c("MYC", "IGF2BP1"),
  c("MYC", "PCBP2"),
  c("MYC", "PKM"),
  c("MYC", "RPL23A"),
  c("MYC", "SERBP1"),
  c("MYC", "SF3B4"),
  c("MYC", "SUB1"),
  c("MYC", "XRN1"),
  c("MYC", "YWHAG")
)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Load expression and metadata files used in the validation analysis.
load_validation_data <- function() {
  list(
    gtex_pheno = fread(GTEX_PHENOTYPE_FILE),
    lihc_surv = fread(LIHC_SURVIVAL_FILE),
    genes_id = fread(GENE_PROBEMAP_FILE),
    liver_expression = readRDS(LIVER_EXPRESSION_FILE)
  )
}

# Extract expression for a gene pair and assign cohort labels.
prepare_pair_expression <- function(expression_matrix, gene1, gene2) {
  gene_indices <- match(c(gene1, gene2), colnames(expression_matrix))
  
  if (any(is.na(gene_indices))) {
    missing_genes <- c(gene1, gene2)[is.na(gene_indices)]
    stop(paste("Missing gene(s) in expression matrix:", paste(missing_genes, collapse = ", ")))
  }
  
  expression_subset <- expression_matrix[, gene_indices, drop = FALSE]
  
  gtex_expression <- expression_subset[grepl("GTEX", rownames(expression_subset)), , drop = FALSE]
  lihc_expression <- expression_subset[grepl("TCGA", rownames(expression_subset)), , drop = FALSE]
  
  colnames(gtex_expression) <- c("gene1", "gene2")
  colnames(lihc_expression) <- c("gene1", "gene2")
  
  gtex_expression <- as.data.frame(gtex_expression)
  lihc_expression <- as.data.frame(lihc_expression)
  
  gtex_expression$group <- "GTEX"
  lihc_expression$group <- "TCGA"
  
  rbind(gtex_expression, lihc_expression)
}

# Create pairwise co-expression plot with regression line and R2/p-value labels.
plot_pair_coexpression <- function(expression_matrix, gene1, gene2) {
  plot_data <- prepare_pair_expression(expression_matrix, gene1, gene2)
  
  ggplot(plot_data, aes(gene1, gene2, colour = group)) +
    geom_point(shape = 16, size = POINT_SIZE, alpha = POINT_ALPHA) +
    stat_poly_line(se = FALSE, linewidth = 2) +
    stat_poly_eq(use_label(c("R2", "p")), size = 5) +
    theme_classic2() +
    scale_color_manual("", values = c("#00AAD8", "#D84500")) +
    theme(text = element_text(size = 16), legend.position = "top") +
    xlab(gene1) +
    ylab(gene2)
}

# Plot and save a grid of co-expression plots.
save_coexpression_grid <- function(expression_matrix, gene_pairs, output_file, ncol, nrow, width, height) {
  plots <- lapply(gene_pairs, function(pair) {
    plot_pair_coexpression(expression_matrix, pair[1], pair[2])
  })
  
  combined_plot <- ggarrange(
    plotlist = plots,
    ncol = ncol,
    nrow = nrow,
    labels = "AUTO",
    common.legend = TRUE
  )
  
  ggsave(
    combined_plot,
    file = output_file,
    limitsize = TRUE,
    width = width,
    height = height,
    dpi = PLOT_DPI,
    bg = "white"
  )
  
  combined_plot
}

# -----------------------------
# LOAD DATA
# -----------------------------
validation_data <- load_validation_data()
liver_expression <- validation_data$liver_expression

# Optional sanity check: number of GTEx and TCGA samples.
n_gtex_samples <- nrow(liver_expression[grepl("GTEX", rownames(liver_expression)), , drop = FALSE])
n_tcga_samples <- nrow(liver_expression[grepl("TCGA", rownames(liver_expression)), , drop = FALSE])

message("GTEx samples: ", n_gtex_samples)
message("TCGA samples: ", n_tcga_samples)

# -----------------------------
# SELECTED PAIR CO-EXPRESSION VALIDATION
# -----------------------------
selected_pair_plot <- save_coexpression_grid(
  expression_matrix = liver_expression,
  gene_pairs = selected_gene_pairs,
  output_file = OUTPUT_COEXPRESSION_FIGURE,
  ncol = 4,
  nrow = 5,
  width = 15,
  height = 20
)

# -----------------------------
# MYC-RELATED CO-EXPRESSION VALIDATION
# -----------------------------
myc_pair_plot <- save_coexpression_grid(
  expression_matrix = liver_expression,
  gene_pairs = myc_gene_pairs,
  output_file = OUTPUT_MYC_COEXPRESSION_FIGURE,
  ncol = 4,
  nrow = 4,
  width = 14,
  height = 18
)

# ============================================================
# END OF SCRIPT
# ============================================================