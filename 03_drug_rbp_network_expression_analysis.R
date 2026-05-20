# ============================================================
# Script: 03_drug_rbp_network_expression_analysis.R
# Purpose: Build drug–RBP interaction networks from CTD/DrugBank
#          records and visualize associated RBP expression changes
#          between LIHC and GTEx liver samples.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
DRUG_DIR <- file.path(INPUT_DIR, "drug_interactions/")
VALIDATION_DIR <- file.path(INPUT_DIR, "validation_UCSC_LIHC_GTEx/")
RESULTS_DIR <- file.path(INPUT_DIR, "results/")
OUTPUT_DIR <- "path/to/output/"

# DrugBank vocabulary and CTD interaction files
DRUGBANK_VOCABULARY_FILE <- file.path(DRUG_DIR, "drugbank_vocabulary.csv")
CTD_FILE_SET_1 <- file.path(DRUG_DIR, "CTD_gene_cgixns_set1.tsv")
CTD_FILE_SET_2 <- file.path(DRUG_DIR, "CTD_gene_cgixns_set2.tsv")

# LIHC/GTEx expression data
GTEX_PHENOTYPE_FILE <- file.path(VALIDATION_DIR, "GTEX_phenotype.gz")
LIHC_SURVIVAL_FILE <- file.path(VALIDATION_DIR, "survival_LIHC_survival.txt")
GENE_PROBEMAP_FILE <- file.path(VALIDATION_DIR, "probeMap_gencode.v23.annotation.gene.probemap")
LIVER_EXPRESSION_FILE <- file.path(RESULTS_DIR, "liver_tcga_gtex.rds")

# Analysis parameters
MIN_DRUG_TARGET_FRACTION <- 0.5
DATASET_1_N_DRUG_NODES <- 6
DATASET_2_N_DRUG_NODES <- 4

# Output files
OUTPUT_FIGURE_SET_1 <- file.path(OUTPUT_DIR, "drug_rbp_network_set1.svg")
OUTPUT_FIGURE_SET_2 <- file.path(OUTPUT_DIR, "drug_rbp_network_set2.svg")

# -----------------------------
# LIBRARIES
# -----------------------------
library(tidyverse)
library(data.table)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggpubr)
library(scales)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Load and filter CTD drug–gene interaction data using DrugBank vocabulary.
load_drug_interactions <- function(ctd_file, drugbank_vocabulary_file, min_target_fraction = 0.5) {
  ctd_table <- read.delim(ctd_file)
  ctd_table <- ctd_table[ctd_table$Organism == "Homo sapiens", ]
  
  drugbank_vocab <- read.csv(drugbank_vocabulary_file)
  ctd_table <- ctd_table[!is.na(match(ctd_table$ChemicalName, drugbank_vocab$Common.name)), ]
  
  # Remove duplicate drug-gene records.
  ctd_table <- ctd_table[!duplicated(paste0(ctd_table$ChemicalName, ctd_table$X..Input)), ]
  
  # Keep drugs linked to a large fraction of genes in the queried interaction set.
  drug_fraction <- sort(table(ctd_table$ChemicalName)) / length(unique(ctd_table$GeneSymbol))
  selected_drugs <- names(drug_fraction[drug_fraction > min_target_fraction])
  
  ctd_table[!is.na(match(ctd_table$ChemicalName, selected_drugs)), c(2, 5, 10)]
}

# Standardize interaction action labels for plotting.
clean_interaction_labels <- function(labels) {
  labels <- gsub("^", " ", labels, fixed = TRUE)
  labels[labels == "affects reaction|increases activity"] <- "increases activity"
  labels[labels == "affects cotreatment|increases expression"] <- "increases expression"
  labels
}

# Build and plot the drug-RBP interaction graph.
plot_drug_rbp_network <- function(drug_interactions, n_drug_nodes) {
  graph_obj <- graph_from_data_frame(drug_interactions, directed = TRUE)
  
  V(graph_obj)$Node_degree <- igraph::degree(graph_obj)
  
  node_type <- rep("RBP", length(V(graph_obj)))
  node_type[seq_len(n_drug_nodes)] <- "drug"
  V(graph_obj)$Node_shape <- node_type
  
  edge_labels <- clean_interaction_labels(drug_interactions$InteractionActions)
  E(graph_obj)$sign <- factor(edge_labels, levels = names(table(edge_labels)))
  
  network_plot <- ggraph(graph_obj, layout = "linear", circular = TRUE) +
    geom_edge_arc(
      width = 1,
      aes(color = sign),
      alpha = 0.6,
      strength = 0.3,
      angle_calc = "along",
      label_dodge = unit(2.5, "mm"),
      arrow = arrow(length = unit(10, "pt")),
      end_cap = circle(20, "pt")
    ) +
    geom_node_point(size = 10, aes(color = Node_shape, shape = Node_shape)) +
    geom_node_text(aes(label = name), colour = "#3F3F3F") +
    scale_color_manual(values = c("#D89FD7", "#D8D8D8"), name = "") +
    scale_edge_color_manual(
      values = c("#BDBDBD", "black", "#FF320D", "#FD8A2F", "#DE41C9", "#34DF4B", "#0D91FF", "#2AEFE3"),
      name = "Interaction"
    ) +
    scale_shape_manual(values = c(18, 19), name = "") +
    theme_graph(base_size = 14) +
    theme(legend.position = "bottom") +
    labs(size = "Node degree", shape = "Node shape")
  
  list(graph = graph_obj, plot = network_plot)
}

# Calculate LIHC versus GTEx liver median expression fold-change for RBP nodes.
calculate_expression_fold_change <- function(graph_obj, liver_expression, n_drug_nodes, flip_fc_indices = integer()) {
  rbp_names <- names(V(graph_obj))[-seq_len(n_drug_nodes)]
  
  # Gene symbol harmonization used in the original analysis.
  rbp_names_for_expression <- rbp_names
  rbp_names_for_expression[rbp_names_for_expression == "RACK1"] <- "GNB2L1"
  
  expression_subset <- liver_expression[, match(rbp_names_for_expression, colnames(liver_expression))]
  expression_subset <- expression_subset[, !is.na(colnames(expression_subset)), drop = FALSE]
  
  gtex_expression <- expression_subset[grepl("GTEX", rownames(expression_subset)), , drop = FALSE]
  lihc_expression <- expression_subset[grepl("TCGA", rownames(expression_subset)), , drop = FALSE]
  
  fold_change_table <- data.frame(
    names = names(log2(apply(lihc_expression, 2, median) / apply(gtex_expression, 2, median))),
    FC = unname(log2(abs(apply(lihc_expression, 2, median)) / abs(apply(gtex_expression, 2, median))))
  )
  
  if (length(flip_fc_indices) > 0) {
    fold_change_table$FC[flip_fc_indices] <- -fold_change_table$FC[flip_fc_indices]
  }
  
  fold_change_table$names[fold_change_table$names == "GNB2L1"] <- "RACK1"
  fold_change_table
}

# Plot RBP expression fold-change.
plot_expression_fold_change <- function(fold_change_table) {
  ggplot(fold_change_table, aes(x = reorder(names, FC), y = FC, fill = FC)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_gradientn(
      colours = c("#32AAD0", "#EBEBEB", "#D03232"),
      values = rescale(c(min(fold_change_table$FC), 0, max(fold_change_table$FC))),
      guide = "colorbar",
      limits = c(min(fold_change_table$FC), max(fold_change_table$FC))
    ) +
    labs(x = "RBP", y = expression("Log"[2] * "FC")) +
    theme_bw() +
    theme(text = element_text(size = 18), legend.position = "none") +
    coord_flip()
}

# Run the complete plotting workflow for one CTD input table.
run_drug_rbp_analysis <- function(ctd_file, n_drug_nodes, output_file, plot_width, expression_width_ratio, flip_fc_indices = integer()) {
  drug_interactions <- load_drug_interactions(
    ctd_file = ctd_file,
    drugbank_vocabulary_file = DRUGBANK_VOCABULARY_FILE,
    min_target_fraction = MIN_DRUG_TARGET_FRACTION
  )
  
  network_result <- plot_drug_rbp_network(
    drug_interactions = drug_interactions,
    n_drug_nodes = n_drug_nodes
  )
  
  fold_change_table <- calculate_expression_fold_change(
    graph_obj = network_result$graph,
    liver_expression = liver_expression,
    n_drug_nodes = n_drug_nodes,
    flip_fc_indices = flip_fc_indices
  )
  
  expression_plot <- plot_expression_fold_change(fold_change_table)
  
  combined_plot <- ggarrange(
    network_result$plot,
    expression_plot,
    ncol = 2,
    labels = "AUTO",
    widths = c(1, expression_width_ratio),
    font.label = list(size = 20)
  )
  
  ggsave(
    filename = output_file,
    plot = combined_plot,
    width = plot_width,
    height = 7
  )
  
  list(
    drug_interactions = drug_interactions,
    network = network_result$graph,
    fold_change = fold_change_table,
    plot = combined_plot
  )
}

# -----------------------------
# LOAD EXPRESSION DATA
# -----------------------------
# These files are kept here for reproducibility, although only the expression
# matrix is used directly in this script.
gtex_pheno <- fread(GTEX_PHENOTYPE_FILE)
lihc_surv <- fread(LIHC_SURVIVAL_FILE)
genes_id <- fread(GENE_PROBEMAP_FILE)
liver_expression <- readRDS(LIVER_EXPRESSION_FILE)

# -----------------------------
# DRUG-RBP NETWORK SET 1
# -----------------------------
result_set_1 <- run_drug_rbp_analysis(
  ctd_file = CTD_FILE_SET_1,
  n_drug_nodes = DATASET_1_N_DRUG_NODES,
  output_file = OUTPUT_FIGURE_SET_1,
  plot_width = 10,
  expression_width_ratio = 0.5,
  flip_fc_indices = 12
)

# -----------------------------
# DRUG-RBP NETWORK SET 2
# -----------------------------
result_set_2 <- run_drug_rbp_analysis(
  ctd_file = CTD_FILE_SET_2,
  n_drug_nodes = DATASET_2_N_DRUG_NODES,
  output_file = OUTPUT_FIGURE_SET_2,
  plot_width = 11,
  expression_width_ratio = 0.6,
  flip_fc_indices = c(9, 14)
)

# ============================================================
# END OF SCRIPT
# ============================================================