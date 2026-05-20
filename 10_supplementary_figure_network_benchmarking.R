# ============================================================
# Script: 10_supplementary_figure_network_benchmarking.R
# Purpose: Generate supplementary figures for benchmarking and
#          comparing network inference methods on simulated
#          ENCODE-like datasets (HepG2/K562).
#
# Supplementary figures include:
#   - Method similarity clustering (Hamming distance)
#   - Activation/inhibition composition heatmaps
#   - Consensus-vs-single method PPV/F1 comparisons
#   - True/false positive loss distributions
#   - Consensus threshold optimization curves
#   - Method contribution pie charts
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/supplementary_figures/"

# Select simulation dataset
CELL_LINE <- "HepG2"   # Options: "HepG2", "K562"

# Input directories
SIMULATION_DIR <- file.path(INPUT_DIR, "ENCODE_simulations/")
NETWORK_DIR <- file.path(SIMULATION_DIR, paste0(tolower(CELL_LINE), "like_nets/"))

# Gold-standard network
GOLD_STANDARD_FILE <- ifelse(
  CELL_LINE == "HepG2",
  file.path(SIMULATION_DIR, "SHEPG2_GRN_SNRL_S4.txt"),
  file.path(SIMULATION_DIR, "SK562_GRN_SNRL_S4.txt")
)

# Output files (supplementary figures)
OUTPUT_CLUSTER_TREE <- file.path(OUTPUT_DIR, paste0("supp_", tolower(CELL_LINE), "_method_clustering.svg"))
OUTPUT_METHOD_SIMILARITY <- file.path(OUTPUT_DIR, paste0("supp_", tolower(CELL_LINE), "_method_similarity_heatmap.svg"))
OUTPUT_CONSENSUS_COMPARISON <- file.path(OUTPUT_DIR, paste0("supp_", tolower(CELL_LINE), "_consensus_comparison.svg"))
OUTPUT_BENCHMARK_SUMMARY <- file.path(OUTPUT_DIR, paste0("supp_", tolower(CELL_LINE), "_benchmarking_summary.svg"))

# Parameters
REMOVE_DIAGONAL <- TRUE
CONSENSUS_CLUSTER_RANGE <- c(3, 9)

# -----------------------------
# LIBRARIES
# -----------------------------
library(R.matlab)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(e1071)
library(ape)
library(pheatmap)
library(ggcorrplot)
library(reshape2)
library(caret)
library(Ckmeans.1d.dp)
library(fmsb)
library(ggplotify)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Convert weighted network into signed adjacency.
binarize_network <- function(net) {
  net[net < 0] <- -1
  net[net > 0] <- 1
  net
}

# Flatten network into binary edge-presence vector.
flatten_binary_network <- function(net) {
  as.numeric(as.logical(as.vector(as.matrix(net))))
}

# Flatten signed network into signed vector.
flatten_signed_network <- function(net) {
  as.vector(as.matrix(net))
}

# Compute pairwise consensus performance.
compute_pairwise_consensus <- function(net1, net2, gold_standard, remove_diagonal = TRUE) {
  if (remove_diagonal) {
    diag_idx <- which(diag(1, nrow(gold_standard)) == 1)
    net1[diag_idx] <- 0
    net2[diag_idx] <- 0
    gold_standard[diag_idx] <- 0
  }
  
  cm1 <- confusionMatrix(as.factor(net1), as.factor(gold_standard))
  cm2 <- confusionMatrix(as.factor(net2), as.factor(gold_standard))
  
  consensus <- as.numeric(as.logical(as.integer(as.factor(net1)) - 1L) &
                            as.logical(as.integer(as.factor(net2)) - 1L))
  cm_consensus <- confusionMatrix(as.factor(consensus), as.factor(gold_standard))
  
  c(
    PPV1 = cm1$table[2,2] / sum(cm1$table[2,]),
    PPV2 = cm2$table[2,2] / sum(cm2$table[2,]),
    PPV_consensus = cm_consensus$table[2,2] / sum(cm_consensus$table[2,]),
    TP1 = cm1$table[2,2],
    TP2 = cm2$table[2,2],
    TP_consensus = cm_consensus$table[2,2],
    FP1 = cm1$table[2,1],
    FP2 = cm2$table[2,1],
    FP_consensus = cm_consensus$table[2,1]
  )
}

# Evaluate a consensus network against the gold standard.
evaluate_network <- function(network, gold_standard) {
  diag_idx <- which(diag(1, nrow(gold_standard)) == 1)
  gold_standard[diag_idx] <- 0
  
  cm <- confusionMatrix(
    as.factor(as.numeric(as.logical(as.vector(gold_standard)))),
    as.factor(as.numeric(as.logical(as.vector(network))))
  )
  
  TP <- cm$table[2,2]
  FP <- cm$table[1,2]
  FN <- cm$table[2,1]
  TN <- cm$table[1,1]
  
  list(
    F1 = 2 * TP / (2 * TP + FP + FN),
    PPV = TP / (TP + FP),
    BACC = ((TP / (TP + FN)) + (TN / (TN + FP))) / 2,
    TP = TP,
    FP = FP
  )
}

# -----------------------------
# LOAD NETWORKS
# -----------------------------
network_files <- list.files(NETWORK_DIR)
method_names <- sapply(strsplit(network_files, "_"), `[`, 2)

gold_standard <- read.table(GOLD_STANDARD_FILE, sep = ",")
gold_standard_signed <- binarize_network(gold_standard)
gold_standard_binary <- flatten_binary_network(gold_standard_signed)

network_list <- lapply(network_files, function(file) {
  binarize_network(read.table(file.path(NETWORK_DIR, file), sep = ","))
})

# -----------------------------
# SUPPLEMENTARY FIGURE 1:
# METHOD CLUSTERING (HAMMING DISTANCE)
# -----------------------------
network_binary_df <- lapply(network_list, flatten_binary_network) %>% bind_cols()
network_binary_df <- cbind(gold_standard_binary, network_binary_df)
colnames(network_binary_df) <- c("GS", method_names)

distance_matrix <- hamming.distance(t(network_binary_df)) %>% as.dist()
hierarchical_clustering <- hclust(distance_matrix, method = "ward.D")

svg(OUTPUT_CLUSTER_TREE, width = 8, height = 8)
plot(as.phylo(hierarchical_clustering))
dev.off()

# -----------------------------
# SUPPLEMENTARY FIGURE 2:
# METHOD SIMILARITY HEATMAP
# -----------------------------
heatmap_matrix <- hamming.distance(t(network_binary_df))

p_heatmap <- pheatmap(
  heatmap_matrix,
  clustering_method = "ward.D",
  fontsize = 14,
  color = colorRampPalette(c("#a561b0", "white", "#58ba54"))(20),
  border_color = "black"
) %>% as.ggplot()

ggsave(OUTPUT_METHOD_SIMILARITY, plot = p_heatmap, width = 8, height = 8)

# -----------------------------
# SUPPLEMENTARY FIGURE 3:
# PAIRWISE CONSENSUS PPV IMPROVEMENT
# -----------------------------
method_pairs <- combn(method_names, 2)
pairwise_results <- list()

ppv_matrix <- matrix(0, nrow = length(method_names), ncol = length(method_names))
rownames(ppv_matrix) <- method_names
colnames(ppv_matrix) <- method_names

for (i in seq_len(ncol(method_pairs))) {
  result <- compute_pairwise_consensus(
    network_binary_df[, method_pairs[1, i]],
    network_binary_df[, method_pairs[2, i]],
    network_binary_df[, "GS"],
    REMOVE_DIAGONAL
  )
  
  pairwise_results[[i]] <- result
  
  ppv_matrix[method_pairs[1, i], method_pairs[2, i]] <- result["PPV_consensus"]
}

p_consensus <- ggcorrplot(
  ppv_matrix,
  method = "circle",
  type = "upper"
) +
  theme(text = element_text(size = 16))

ggsave(OUTPUT_CONSENSUS_COMPARISON, plot = p_consensus, width = 8, height = 8)

# -----------------------------
# SUPPLEMENTARY FIGURE 4:
# CONSENSUS THRESHOLD OPTIMIZATION
# -----------------------------
absolute_networks <- lapply(network_list, abs)
average_network <- Reduce("+", absolute_networks) / length(absolute_networks)
thresholds <- sort(unique(abs(average_network[average_network != 0])))

consensus_networks <- lapply(thresholds, function(threshold) {
  net <- average_network
  net[abs(net) < threshold] <- 0
  net[net != 0] <- 1
  net
})

benchmark_results <- lapply(consensus_networks, evaluate_network, gold_standard = gold_standard_signed)

benchmark_df <- data.frame(
  threshold = seq_along(thresholds),
  F1 = sapply(benchmark_results, `[[`, "F1"),
  PPV = sapply(benchmark_results, `[[`, "PPV")
)

p_f1 <- ggplot(benchmark_df, aes(x = threshold, y = F1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  ylab("F1 score") +
  xlab("Consensus threshold")

p_ppv <- ggplot(benchmark_df, aes(x = threshold, y = PPV)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  ylab("PPV") +
  xlab("Consensus threshold")

benchmark_summary <- ggarrange(
  p_f1,
  p_ppv,
  labels = c("A", "B"),
  ncol = 2
)

ggsave(OUTPUT_BENCHMARK_SUMMARY, plot = benchmark_summary, width = 12, height = 6)

# ============================================================
# END OF SCRIPT
# ============================================================