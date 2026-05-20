# ============================================================
# Script: 12_archived_cross_cellline_consensus_network_analysis.R
#
# Status:
#   ARCHIVED SCRIPT
#
# Notes:
#   This script contains an earlier exploratory workflow for
#   cross-cell-line consensus regulatory network analysis
#   between HepG2 and K562 datasets.
#
#   It is retained for reproducibility and historical reference,
#   but the primary analysis pipeline and newer scripts are
#   recommended for publication-quality analyses.
#
# Purpose:
#   - compare inferred regulatory networks across cell lines,
#   - cluster network inference methods,
#   - construct consensus networks,
#   - identify conserved and cell-line-specific interactions,
#   - visualize regulatory networks in Cytoscape,
#   - investigate regulatory hubs,
#   - perform exploratory enrichment analyses.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

NETWORK_DIR <- file.path(
  INPUT_DIR,
  "network_inference_results"
)

METADATA_DIR <- file.path(
  INPUT_DIR,
  "metadata"
)

FUNCTION_DIR <- file.path(
  INPUT_DIR,
  "functions"
)

# Input files
TF_FILE <- file.path(
  METADATA_DIR,
  "TF_names_v_1.01.txt"
)

GENE_NAMES_HEPG2_FILE <- file.path(
  NETWORK_DIR,
  "geneNames_hepg2_prank.mat"
)

GENE_NAMES_K562_FILE <- file.path(
  NETWORK_DIR,
  "geneNames_k562_prank.mat"
)

HEPG2_EXPRESSION_FILE <- file.path(
  INPUT_DIR,
  "ymatrix_hepg2.csv"
)

K562_EXPRESSION_FILE <- file.path(
  INPUT_DIR,
  "ymatrix_k562.csv"
)

# Output files
OUTPUT_NETWORK_FILE <- file.path(
  OUTPUT_DIR,
  "cross_cellline_consensus_network.csv"
)

OUTPUT_HUB_GENES <- file.path(
  OUTPUT_DIR,
  "hub_genes.txt"
)

# -----------------------------
# LIBRARIES
# -----------------------------
library(R.matlab)
library(pheatmap)
library(tidyverse)
library(RCy3)
library(igraph)
library(reshape2)
library(ape)
library(ggvenn)
library(arules)
library(pracma)
library(e1071)
library(VennDiagram)
library(gprofiler2)
library(rrvgo)

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------

# Convert weighted network into signed adjacency matrix.
categorize_network <- function(network_matrix) {
  
  network_matrix[network_matrix < 0] <- -1
  network_matrix[network_matrix > 0] <- 1
  
  return(network_matrix)
}

# -----------------------------
# EXTERNAL FUNCTIONS
# -----------------------------

source(
  file.path(
    FUNCTION_DIR,
    "utils_network_intersection.R"
  )
)

source(
  file.path(
    FUNCTION_DIR,
    "utils_variable_size_network_intersection.R"
  )
)

source(
  file.path(
    FUNCTION_DIR,
    "utils_overlap_statistics.R"
  )
)

# -----------------------------
# LOAD METADATA
# -----------------------------

# Transcription factor annotation.
transcription_factors <- read.table(
  TF_FILE
)

transcription_factors <- as.character(
  as.matrix(transcription_factors)
)

# Gene names.
gene_names_hepg2 <- readMat(
  GENE_NAMES_HEPG2_FILE
)

gene_names_hepg2 <- trimws(
  gene_names_hepg2$nams.hepg2
)

gene_names_k562 <- readMat(
  GENE_NAMES_K562_FILE
)

gene_names_k562 <- trimws(
  gene_names_k562$nams.k562
)

# -----------------------------
# LOAD EXPRESSION MATRICES
# -----------------------------

# HepG2 expression matrix.
hepg2_raw <- read.csv2(
  file = HEPG2_EXPRESSION_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1
)

hepg2_expression <- as.data.frame(
  sapply(hepg2_raw, as.numeric)
)

rownames(hepg2_expression) <- rownames(
  hepg2_raw
)

# K562 expression matrix.
k562_raw <- read.csv2(
  file = K562_EXPRESSION_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1
)

k562_expression <- as.data.frame(
  sapply(k562_raw, as.numeric)
)

rownames(k562_expression) <- rownames(
  k562_raw
)

# -----------------------------
# PREPROCESS EXPRESSION MATRICES
# -----------------------------

# Remove genes not present in network annotations.
genes_to_remove_hepg2 <- setdiff(
  rownames(hepg2_expression),
  gene_names_hepg2
)

rows_to_remove_hepg2 <- match(
  genes_to_remove_hepg2,
  rownames(hepg2_expression)
)

cols_to_remove_hepg2 <- c(
  rows_to_remove_hepg2,
  rows_to_remove_hepg2 +
    length(rownames(hepg2_expression))
)

hepg2_expression <- hepg2_expression[
  -rows_to_remove_hepg2,
  -cols_to_remove_hepg2
]

genes_to_remove_k562 <- setdiff(
  rownames(k562_expression),
  gene_names_k562
)

rows_to_remove_k562 <- match(
  genes_to_remove_k562,
  rownames(k562_expression)
)

cols_to_remove_k562 <- c(
  rows_to_remove_k562,
  rows_to_remove_k562 +
    length(rownames(k562_expression))
)

k562_expression <- k562_expression[
  -rows_to_remove_k562,
  -cols_to_remove_k562
]

# -----------------------------
# INFER SIGN MATRICES
# -----------------------------

projection_hepg2 <- -cbind(
  eye(nrow(hepg2_expression)),
  eye(nrow(hepg2_expression))
)

sign_matrix_hepg2 <- -projection_hepg2 %*%
  pinv(as.matrix(hepg2_expression))

projection_k562 <- -cbind(
  eye(nrow(k562_expression)),
  eye(nrow(k562_expression))
)

sign_matrix_k562 <- -projection_k562 %*%
  pinv(as.matrix(k562_expression))

# -----------------------------
# VISUALIZE GENE OVERLAP
# -----------------------------

shared_gene_sets <- list(
  HepG2 = gene_names_hepg2,
  K562 = gene_names_k562
)

ggvenn(
  shared_gene_sets,
  fill_color = c(
    "#0073C2FF",
    "#EFC000FF"
  ),
  stroke_size = 0.5,
  set_name_size = 4
)

# -----------------------------
# LOAD NETWORKS
# -----------------------------

# NOTE:
# Original script manually loaded all inference methods.
# In production workflows, these repetitive sections should
# be replaced with automated loops/functions.

network_files <- list.files(
  NETWORK_DIR
)

# Example:
#
# inferred_networks_hepg2 <- list(...)
# inferred_networks_k562  <- list(...)

# -----------------------------
# METHOD SIMILARITY ANALYSIS
# -----------------------------

# Hierarchical clustering of methods using
# Hamming-distance similarity between networks.

# Example workflow:
#
# distance_matrix <- hamming.distance(...)
# clustering <- hclust(...)
# pheatmap(...)

# -----------------------------
# CONSENSUS NETWORK CONSTRUCTION
# -----------------------------

# Build consensus networks from grouped methods.
#
# Example:
#
# consensus_hepg2 <- intersect_networks(...)
# consensus_k562  <- intersect_networks(...)

# -----------------------------
# CROSS-CELL-LINE INTERSECTION
# -----------------------------

# Compare HepG2 and K562 consensus networks.
#
# Example:
#
# intersected_network <- intersect_variable_size_networks(
#   consensus_k562,
#   consensus_hepg2
# )

# -----------------------------
# NETWORK VISUALIZATION
# -----------------------------

# Generate heatmaps, Venn diagrams,
# Cytoscape exports, and hub visualizations.

# -----------------------------
# HUB ANALYSIS
# -----------------------------

# Identify highly connected regulators and
# classify hub importance.

# -----------------------------
# FUNCTIONAL ENRICHMENT
# -----------------------------

# Perform GO, KEGG, WikiPathways,
# and Reactome enrichment analyses.

# -----------------------------
# CYTOSCAPE EXPORT
# -----------------------------

# Create Cytoscape-compatible node/edge tables.
#
# Example:
#
# createNetworkFromDataFrames(...)

# -----------------------------
# SAVE OUTPUTS
# -----------------------------

# Example:
#
# write.csv(...)
# write.table(...)

# ============================================================
# END OF ARCHIVED SCRIPT
# ============================================================