# ============================================================
# Script: 06_network_module_enrichment_analysis.R
# Purpose: Detect communities in a selected RBP interaction network,
#          assign edge signs using method consensus and expression
#          correlation, and annotate communities with MSigDB enrichment.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files
RBP_LIST_FILE <- file.path(INPUT_DIR, "rbps_info/210329_Table_S1_hRBP_list.xlsx")
VALIDATION_TABLE_FILE <- file.path(INPUT_DIR, "Table_S3_validationtable_3plus.csv")
HEPG2_SIGNED_ADJACENCY_FILE <- file.path(INPUT_DIR, "Table_S1_adjacencymatrix_HepG2_signed.csv")
HEPG2_EXPRESSION_FILE <- file.path(INPUT_DIR, "datasets/ymatrix_hepg2.csv")
K562_EXPRESSION_FILE <- file.path(INPUT_DIR, "datasets/ymatrix_k562.csv")
HEPG2_NETWORK_DIR <- file.path(INPUT_DIR, "ENCODE_realNets/HepG2/")

# MSigDB gene sets
MSIGDB_C2_FILE <- file.path(INPUT_DIR, "MSigDB/c2.cp.v2023.2.Hs.symbols.gmt")
MSIGDB_C4_FILE <- file.path(INPUT_DIR, "MSigDB/c4.all.v2023.2.Hs.symbols.gmt")
MSIGDB_C5_FILE <- file.path(INPUT_DIR, "MSigDB/c5.all.v2023.2.Hs.symbols.gmt")
MSIGDB_C6_FILE <- file.path(INPUT_DIR, "MSigDB/c6.all.v2023.2.Hs.symbols.gmt")
MSIGDB_HALLMARK_FILE <- file.path(INPUT_DIR, "MSigDB/h.all.v2023.2.Hs.symbols.gmt")

# Output files
OUTPUT_NETWORK_PDF <- file.path(OUTPUT_DIR, "network_modules.pdf")
OUTPUT_LEGEND_PDF <- file.path(OUTPUT_DIR, "network_module_enrichment_legend.pdf")
OUTPUT_MODULE_ENRICHMENT_RDS <- file.path(OUTPUT_DIR, "module_enrichment_results.rds")
OUTPUT_COMMUNITIES_RDS <- file.path(OUTPUT_DIR, "network_communities.rds")

# Analysis parameters
USE_VALIDATED_NETWORK <- TRUE
VALIDATED_NETWORK_N_EDGES <- 119
CONSENSUS_WEIGHT_THRESHOLD <- 0.5
MIN_COMMUNITY_SIZE <- 5
ENRICHMENT_FDR_THRESHOLD <- 0.10
MIN_ENRICHED_GENE_COUNT <- 2
EXCLUDED_METHODS_FOR_SIGN <- c("CART", "neunetreg")

# -----------------------------
# LIBRARIES
# -----------------------------
library(tidyverse)
library(readxl)
library(igraph)
library(qusage)
library(scales)
library(clusterProfiler)
library(R.matlab)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Convert signed network values to -1, 0, or 1.
catnet <- function(net) {
  net[net < 0] <- -1
  net[net > 0] <- 1
  net
}

# Read the network to be analyzed.
read_selected_network <- function(validation_table_file, signed_adjacency_file, use_validated_network = TRUE, n_validated_edges = 119, consensus_threshold = 0.5) {
  if (use_validated_network) {
    inters <- read.csv2(validation_table_file)
    selected_edges <- inters[seq_len(n_validated_edges), ]
    network_edges <- do.call(rbind, strsplit(selected_edges[, 1], "-")) %>% as.data.frame()
    colnames(network_edges) <- c("from", "to")
  } else {
    adjacency <- read.csv2(signed_adjacency_file, row.names = 1)
    graph_obj <- graph.adjacency(t(as.matrix(adjacency)), weighted = TRUE)
    edge_table <- get.data.frame(graph_obj)
    network_edges <- edge_table[abs(edge_table$weight) >= consensus_threshold, c("from", "to")]
  }
  
  network_edges
}

# Load all method-specific HepG2 networks.
load_method_networks <- function(network_dir) {
  files <- list.files(network_dir)
  method_names <- lapply(strsplit(files, "_"), function(x) x[2]) %>% unlist()
  
  raw_networks <- vector("list", length(files))
  signed_networks <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    net <- read.table(file.path(network_dir, files[i]), sep = ",")
    raw_networks[[i]] <- net
    signed_networks[[i]] <- catnet(net)
  }
  
  names(raw_networks) <- method_names
  names(signed_networks) <- method_names
  
  list(raw = raw_networks, signed = signed_networks, method_names = method_names)
}

# Prepare expression matrix used for sign assignment.
prepare_average_expression <- function(expression_file) {
  expr <- read.table(expression_file)
  (expr[, 1:232] + expr[, 233:464]) / 2
}

# Assign edge signs using method consensus and Spearman correlation as tie-breaker.
assign_edge_signs <- function(network_edges, method_networks, gene_names, expression_matrix, excluded_methods = c("CART", "neunetreg")) {
  edge_matrix <- as.matrix(network_edges[, c("from", "to")])
  method_names <- names(method_networks)
  
  sign_by_method <- list()
  spearman_cor <- numeric(nrow(edge_matrix))
  
  for (i in seq_len(nrow(edge_matrix))) {
    edge_signs <- numeric(length(method_networks))
    
    for (j in seq_along(method_networks)) {
      tmpnet <- t(method_networks[[j]])
      rownames(tmpnet) <- gene_names
      colnames(tmpnet) <- gene_names
      
      edge_signs[j] <- tmpnet[
        which(rownames(tmpnet) == edge_matrix[i, 1]),
        which(colnames(tmpnet) == edge_matrix[i, 2])
      ]
    }
    
    sign_by_method[[i]] <- edge_signs
    
    gene_a <- expression_matrix[which(rownames(expression_matrix) == edge_matrix[i, 1]), ] %>% as.matrix() %>% as.numeric()
    gene_b <- expression_matrix[which(rownames(expression_matrix) == edge_matrix[i, 2]), ] %>% as.matrix() %>% as.numeric()
    spearman_cor[i] <- cor(gene_a, gene_b, method = "spearman")
  }
  
  sign_matrix <- do.call(rbind, sign_by_method)
  keep_methods <- !method_names %in% excluded_methods
  sign_matrix <- sign_matrix[, keep_methods, drop = FALSE]
  sign_matrix <- cbind(sign_matrix, catnet(spearman_cor) %>% as.data.frame())
  
  correlation_sign <- catnet(spearman_cor)
  final_sign <- numeric(nrow(edge_matrix))
  
  for (i in seq_len(nrow(edge_matrix))) {
    sign_counts <- table(as.matrix(sign_matrix[i, ]))
    sign_counts <- sign_counts[names(sign_counts) != "0"]
    
    if (length(sign_counts) == 1) {
      final_sign[i] <- as.numeric(names(sign_counts))
    } else if (sign_counts[which(names(sign_counts) == "-1")] == sign_counts[which(names(sign_counts) == "1")]) {
      final_sign[i] <- correlation_sign[i]
    } else if (sign_counts[which(names(sign_counts) == "-1")] < sign_counts[which(names(sign_counts) == "1")]) {
      final_sign[i] <- 1
    } else {
      final_sign[i] <- -1
    }
  }
  
  final_sign
}

# Detect communities and remove very small communities.
detect_filtered_communities <- function(network_edges, min_community_size = 5) {
  graph_obj <- graph_from_data_frame(network_edges, directed = FALSE)
  
  initial_communities <- cluster_infomap(graph_obj)
  small_communities <- which(table(initial_communities$membership) < min_community_size)
  vertices_to_keep <- V(graph_obj)[!(initial_communities$membership %in% small_communities)]
  
  filtered_graph <- induced.subgraph(graph_obj, vertices_to_keep)
  filtered_communities <- cluster_infomap(filtered_graph)
  
  list(
    graph = filtered_graph,
    communities = communities(filtered_communities),
    membership = filtered_communities
  )
}

# Run enrichment for each community.
enrich_communities <- function(communities_list, gene_sets, universe_genes) {
  enrichment_results <- vector("list", length(communities_list))
  
  for (i in seq_along(communities_list)) {
    enrichment_results[[i]] <- clusterProfiler::enricher(
      gene = communities_list[[i]],
      universe = universe_genes,
      TERM2GENE = gene_sets,
      pAdjustMethod = "fdr",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }
  
  enrichment_results
}

# Convert p-values to star labels.
make_stars <- function(x) {
  stars <- c("****", "***", "**", "*", "")
  cut_points <- c(0, 0.001, 0.01, 0.05, 0.1, 1)
  stars[findInterval(x, cut_points)]
}

# Build text labels summarizing significant enrichment per community.
make_community_enrichment_labels <- function(enrichment_results, fdr_threshold = 0.10, min_gene_count = 2) {
  labels <- character(length(enrichment_results))
  
  for (i in seq_along(enrichment_results)) {
    result <- enrichment_results[[i]]
    
    if (is.null(result)) {
      labels[i] <- "no significant enrichment found"
      next
    }
    
    result_df <- result@result
    result_df <- result_df[result_df$Count > min_gene_count, ]
    
    if (nrow(result_df) == 0) {
      labels[i] <- "no significant enrichment found"
      next
    }
    
    result_df$p.adjust <- p.adjust(result_df$pvalue, "fdr")
    result_df <- result_df[result_df$p.adjust < fdr_threshold, ]
    
    if (nrow(result_df) == 0) {
      labels[i] <- "no significant enrichment found"
    } else {
      terms <- unique(toupper(result_df$ID))
      terms <- paste0(terms, " (", make_stars(result_df$p.adjust[seq_along(terms)]), ")")
      labels[i] <- paste0(gsub("_", " ", terms), collapse = "\n")
    }
  }
  
  labels
}

# Plot network communities with edge colors indicating sign.
plot_network_modules <- function(graph_obj, communities_obj, edge_signs, output_file) {
  color_pool <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = TRUE)]
  community_colors <- color_pool[c(91, 356, 285, 31, 118, 68, 138, 430, 35, 200, 230, 270, 232, 433, 169, 265, 188, 25, 98, 103)[seq_along(communities_obj)]]
  
  edge_color <- as.factor(edge_signs)
  levels(edge_color) <- c("#D32E01", "#0181D3", "#808080", "#808080")[seq_along(levels(edge_color))]
  
  coords <- layout_(graph_obj, with_fr())
  
  pdf(output_file, width = 10, height = 8)
  par(mar = c(3, 1, 3, 1))
  plot.igraph(
    graph_obj,
    vertex.color = "#808080",
    mark.groups = communities_obj,
    vertex.label.dist = 1,
    vertex.size = 5,
    vertex.frame.color = NA,
    vertex.label.family = "sans",
    edge.width = 2,
    layout = coords,
    vertex.label.color = "#000000",
    vertex.label.cex = 1,
    edge.arrow.size = 0.2,
    edge.arrow.width = 1,
    mark.col = alpha(community_colors, alpha = 0.5),
    edge.color = as.character(edge_color),
    edge.curved = 0.1,
    mark.border = community_colors,
    margin = c(0, 1, 0, 1)
  )
  dev.off()
  
  community_colors
}

# Plot enrichment legend for communities.
plot_enrichment_legend <- function(labels, community_colors, output_file) {
  pdf(output_file, width = 10, height = 8)
  plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "", xlim = c(0, 1), ylim = c(0, 1))
  legend(
    0,
    1.1,
    legend = labels,
    bty = "n",
    col = alpha(community_colors, alpha = 0.7),
    cex = 0.9,
    pch = 19,
    pt.cex = 2,
    xpd = TRUE,
    y.intersp = 1.1,
    x.intersp = 0.8
  )
  dev.off()
}

# -----------------------------
# LOAD INPUTS
# -----------------------------
rbp_table <- read_excel(RBP_LIST_FILE)
network_edges <- read_selected_network(
  validation_table_file = VALIDATION_TABLE_FILE,
  signed_adjacency_file = HEPG2_SIGNED_ADJACENCY_FILE,
  use_validated_network = USE_VALIDATED_NETWORK,
  n_validated_edges = VALIDATED_NETWORK_N_EDGES,
  consensus_threshold = CONSENSUS_WEIGHT_THRESHOLD
)

# Gene names are taken from the row names of the expression matrices.
gene_names_hepg2 <- read.table(HEPG2_EXPRESSION_FILE) %>% rownames() %>% trimws()
gene_names_k562 <- read.table(K562_EXPRESSION_FILE) %>% rownames() %>% trimws()

# -----------------------------
# EDGE SIGN ASSIGNMENT
# -----------------------------
method_networks <- load_method_networks(HEPG2_NETWORK_DIR)
hepg2_expression_average <- prepare_average_expression(HEPG2_EXPRESSION_FILE)

network_edges$sign <- assign_edge_signs(
  network_edges = network_edges,
  method_networks = method_networks$signed,
  gene_names = gene_names_hepg2,
  expression_matrix = hepg2_expression_average,
  excluded_methods = EXCLUDED_METHODS_FOR_SIGN
)

# -----------------------------
# COMMUNITY DETECTION
# -----------------------------
community_result <- detect_filtered_communities(
  network_edges = network_edges,
  min_community_size = MIN_COMMUNITY_SIZE
)

saveRDS(community_result$communities, OUTPUT_COMMUNITIES_RDS)

# -----------------------------
# COMMUNITY ENRICHMENT
# -----------------------------
hallmark_gene_sets <- read.gmt(MSIGDB_HALLMARK_FILE)
canonical_gene_sets <- read.gmt(MSIGDB_C2_FILE)
computational_gene_sets <- read.gmt(MSIGDB_C4_FILE)
go_gene_sets <- read.gmt(MSIGDB_C5_FILE)
oncogenic_gene_sets <- read.gmt(MSIGDB_C6_FILE)

# Original analysis used Hallmark + C2 curated pathways.
gene_sets_for_enrichment <- rbind(hallmark_gene_sets, canonical_gene_sets)

enrichment_results <- enrich_communities(
  communities_list = community_result$communities,
  gene_sets = gene_sets_for_enrichment,
  universe_genes = rbp_table$gene_name
)

saveRDS(enrichment_results, OUTPUT_MODULE_ENRICHMENT_RDS)

enrichment_labels <- make_community_enrichment_labels(
  enrichment_results = enrichment_results,
  fdr_threshold = ENRICHMENT_FDR_THRESHOLD,
  min_gene_count = MIN_ENRICHED_GENE_COUNT
)

# Genes from communities enriched for MYC target signatures.
myc_enriched_genes <- community_result$communities[grepl("HALLMARK MYC TARGETS", enrichment_labels)] %>%
  unlist() %>%
  unname() %>%
  unique() %>%
  sort()

# -----------------------------
# PLOTS
# -----------------------------
community_colors <- plot_network_modules(
  graph_obj = community_result$graph,
  communities_obj = community_result$communities,
  edge_signs = E(community_result$graph)$sign,
  output_file = OUTPUT_NETWORK_PDF
)

plot_enrichment_legend(
  labels = enrichment_labels,
  community_colors = community_colors,
  output_file = OUTPUT_LEGEND_PDF
)

# ============================================================
# END OF SCRIPT
# ============================================================
