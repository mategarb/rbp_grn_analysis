# ============================================================
# Script: 07_network_topology_comparison.R
# Purpose: Compare topology properties of inferred ENCODE networks
#          across cell lines, inference methods, and sparsity criteria,
#          using TRRUST as an external reference network.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input directories and files
NETWORK_DIR <- file.path(INPUT_DIR, "optimal_inferred_networks_encode/")
TRRUST_FILE <- file.path(INPUT_DIR, "reference_networks/trrust_rawdata_human.tsv")

# Output files
OUTPUT_TOPOLOGY_TABLE <- file.path(OUTPUT_DIR, "network_topology_metrics.csv")
OUTPUT_TOPOLOGY_FIGURE <- file.path(OUTPUT_DIR, "network_topology_comparison.svg")

# Network metadata
CELL_LINES <- c("HepG2", "K562")
METHODS <- c("GENIE3", "lasso", "LSCON", "Zscore")
SPARSITY_METRICS <- c("corr", "Q")
SPARSITY_LABELS <- c(corr = "LL", Q = "GOF")
NETWORK_SIZE_LABEL <- "1000"

# Plot parameters
FIGURE_WIDTH <- 8
FIGURE_HEIGHT <- 8
FIGURE_DPI <- 300

# -----------------------------
# LIBRARIES
# -----------------------------
library(tidyverse)
library(igraph)
library(R.matlab)
library(ggplot2)
library(ggpubr)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Build the expected .mat file path for a network.
get_network_file <- function(method, cell_line, sparsity_metric) {
  file_name <- paste0(
    "A_", method, "_", NETWORK_SIZE_LABEL, "_", cell_line, "_", sparsity_metric, ".mat"
  )
  file.path(NETWORK_DIR, file_name)
}

# Read an adjacency matrix from a MATLAB .mat file.
read_network_matrix <- function(file_path, sparsity_metric) {
  mat_object <- readMat(file_path)
  
  matrix_name <- ifelse(sparsity_metric == "corr", "A.corr", "A.Q")
  
  if (!matrix_name %in% names(mat_object)) {
    stop(paste("Expected matrix", matrix_name, "not found in", file_path))
  }
  
  mat_object[[matrix_name]]
}

# Convert network to binary adjacency and calculate topology metrics.
calculate_network_properties <- function(network_matrix) {
  binary_network <- network_matrix
  binary_network[binary_network != 0] <- 1
  
  graph_obj <- graph_from_adjacency_matrix(
    as.matrix(binary_network),
    mode = "directed",
    diag = FALSE
  )
  
  outdegree <- igraph::degree(graph_obj, mode = "out")
  degree_table <- table(outdegree)
  
  log_degree <- log(as.numeric(names(degree_table)) + 1)
  log_frequency <- log(as.numeric(degree_table) + 1)
  
  # R-squared from log-log outdegree distribution.
  degree_fit <- summary(lm(log_degree ~ log_frequency))
  
  data.frame(
    R2 = degree_fit$r.squared,
    mean_outdegree = mean(outdegree),
    betweenness_variance = var(betweenness(graph_obj))
  )
}

# Calculate metrics for all inferred networks.
calculate_all_network_metrics <- function() {
  results <- list()
  
  counter <- 1
  for (cell_line in CELL_LINES) {
    for (sparsity_metric in SPARSITY_METRICS) {
      for (method in METHODS) {
        network_file <- get_network_file(method, cell_line, sparsity_metric)
        network_matrix <- read_network_matrix(network_file, sparsity_metric)
        metrics <- calculate_network_properties(network_matrix)
        
        metrics$method <- method
        metrics$sparsity_metric <- SPARSITY_LABELS[[sparsity_metric]]
        metrics$cell_line <- cell_line
        metrics$file <- basename(network_file)
        
        results[[counter]] <- metrics
        counter <- counter + 1
      }
    }
  }
  
  bind_rows(results)
}

# Load TRRUST as a reference regulatory network and calculate topology metrics.
calculate_trrust_reference_metrics <- function(trrust_file) {
  trrust <- read.table(trrust_file, sep = "\t", header = FALSE)
  
  trrust_edges <- as.matrix(trrust[, 1:2])
  trrust_graph <- graph_from_edgelist(trrust_edges, directed = TRUE)
  trrust_matrix <- as_adjacency_matrix(trrust_graph, sparse = FALSE)
  
  calculate_network_properties(trrust_matrix)
}

# Plot R-squared of log-log degree distribution.
plot_degree_distribution_r2 <- function(metrics_df, reference_metrics) {
  ggplot(metrics_df, aes(fill = sparsity_metric, y = R2, x = method)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_hline(yintercept = reference_metrics$R2, linetype = "dashed", color = "#424242", linewidth = 1) +
    facet_wrap(~cell_line) +
    theme_classic() +
    theme(text = element_text(size = 18)) +
    scale_fill_manual(values = c("#E1992A", "#13A79C"), name = "sparsity metric") +
    xlab("") +
    ylab("R-squared of\nlog-log distribution")
}

# Plot mean outdegree.
plot_mean_outdegree <- function(metrics_df, reference_metrics) {
  ggplot(metrics_df, aes(fill = sparsity_metric, y = mean_outdegree, x = method)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_hline(yintercept = reference_metrics$mean_outdegree, linetype = "dashed", color = "#424242", linewidth = 1) +
    facet_wrap(~cell_line) +
    theme_classic() +
    theme(text = element_text(size = 18)) +
    scale_fill_manual(values = c("#E1992A", "#13A79C"), name = "sparsity metric") +
    xlab("") +
    ylab("Mean\noutdegree")
}

# Optional plot: variance in betweenness centrality.
plot_betweenness_variance <- function(metrics_df, reference_metrics) {
  ggplot(metrics_df, aes(fill = sparsity_metric, y = betweenness_variance, x = method)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_hline(yintercept = reference_metrics$betweenness_variance, linetype = "dashed", color = "#424242", linewidth = 1) +
    facet_wrap(~cell_line) +
    theme_classic() +
    theme(text = element_text(size = 18)) +
    scale_fill_manual(values = c("#E1992A", "#13A79C"), name = "sparsity metric") +
    xlab("") +
    scale_y_log10() +
    ylab("Log betweenness\ncentrality variance")
}

# -----------------------------
# CALCULATE NETWORK METRICS
# -----------------------------
topology_metrics <- calculate_all_network_metrics()
trrust_reference_metrics <- calculate_trrust_reference_metrics(TRRUST_FILE)

write.csv(topology_metrics, OUTPUT_TOPOLOGY_TABLE, row.names = FALSE)

# -----------------------------
# PLOT TOPOLOGY COMPARISON
# -----------------------------
p_r2 <- plot_degree_distribution_r2(topology_metrics, trrust_reference_metrics)
p_outdegree <- plot_mean_outdegree(topology_metrics, trrust_reference_metrics)

combined_plot <- ggarrange(
  p_r2,
  p_outdegree,
  labels = "AUTO",
  common.legend = TRUE,
  nrow = 2,
  align = "v"
)

ggsave(
  filename = OUTPUT_TOPOLOGY_FIGURE,
  plot = combined_plot,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = FIGURE_DPI,
  bg = "white"
)

# ============================================================
# END OF SCRIPT
# ============================================================
