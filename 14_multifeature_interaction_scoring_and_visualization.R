# ============================================================
# Script: 14_multifeature_interaction_scoring_and_visualization.R
# Purpose: Integrate multiple validation features to score,
#          rank, classify, and visualize regulatory interactions.
#
# This script:
#   - trains a classification model for interaction confidence,
#   - predicts confidence for unlabeled interactions,
#   - integrates multiple biological evidence layers,
#   - computes interaction ranking scores,
#   - visualizes top-ranked interactions,
#   - generates interaction network figures,
#   - performs disease enrichment analyses.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files
EDGE_FEATURE_FILE <- file.path(
  INPUT_DIR,
  "inters_hepg2_eclip.rds"
)

RBP_ANNOTATION_FILE <- file.path(
  INPUT_DIR,
  "rbps_info",
  "210329_Table_S1_hRBP_list.xlsx"
)

# Output files
OUTPUT_FINAL_EDGES <- file.path(
  OUTPUT_DIR,
  "edges_final_hepg2.rds"
)

OUTPUT_RANKED_NETWORK <- file.path(
  OUTPUT_DIR,
  "net_scores_val.rds"
)

OUTPUT_VALIDATION_HEATMAP <- file.path(
  OUTPUT_DIR,
  "top_interaction_validation_heatmap.svg"
)

OUTPUT_NETWORK_VISUALIZATION <- file.path(
  OUTPUT_DIR,
  "top_interaction_network.svg"
)

OUTPUT_COMBINED_FIGURE <- file.path(
  OUTPUT_DIR,
  "interaction_summary_figure.svg"
)

OUTPUT_TOP_NETWORK <- file.path(
  OUTPUT_DIR,
  "hepg2_regulome_top64.csv"
)

# Parameters
TOP_INTERACTIONS <- 64
N_CLUSTERS <- 3

# -----------------------------
# LIBRARIES
# -----------------------------
library(classInt)
library(caTools)
library(party)
library(dplyr)
library(magrittr)
library(ISLR)
library(rpart)
library(rpart.plot)
library(tidyr)
library(tidyverse)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(ggpubr)
library(gtable)
library(ggnetwork)
library(ggnet)
library(scales)
library(mclust)
library(caret)
library(igraph)
library(gprofiler2)
library(readxl)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Cluster non-zero values while preserving NA/zero values.
cluster_with_missing <- function(values, clusters) {
  
  values_copy <- values
  
  non_missing <- which(!is.na(values_copy))
  zero_indices <- which(values_copy == 0)
  
  output_vector <- rep(
    NA,
    length(values_copy)
  )
  
  clustering <- Mclust(
    as.numeric(
      na.omit(values_copy[-zero_indices])
    ),
    clusters
  )$classification
  
  output_vector[
    setdiff(non_missing, zero_indices)
  ] <- clustering
  
  output_vector[zero_indices] <- 0
  
  return(output_vector)
}

# Convert significance levels into symbols.
convert_significance_labels <- function(values) {
  
  output <- values
  
  for (i in seq_along(values)) {
    
    if (is.na(values[i])) {
      
      output[i] <- "NA"
      
    } else if (values[i] == 0) {
      
      output[i] <- "ns"
      
    } else if (values[i] == 1) {
      
      output[i] <- "*"
      
    } else if (values[i] == 2) {
      
      output[i] <- "**"
      
    } else if (values[i] == 3) {
      
      output[i] <- "***"
      
    } else if (values[i] == 4) {
      
      output[i] <- "****"
    }
  }
  
  return(output)
}

# -----------------------------
# LOAD FEATURE TABLE
# -----------------------------

edge_features <- readRDS(
  EDGE_FEATURE_FILE
)

edge_features$cfreq <- abs(
  edge_features$cfreq
)

# -----------------------------
# PREPARE TRAINING DATA
# -----------------------------

training_edges <- edge_features[
  !is.na(edge_features$fc_score),
]

training_scores <- training_edges$fc_score

training_edges <- training_edges[
  ,
  -match(
    c(
      "fc_score",
      "fc_score_binned"
    ),
    colnames(training_edges)
  )
]

binary_scores <- training_scores

binary_scores[
  training_scores < mean(training_scores)
] <- "low"

binary_scores[
  training_scores >= mean(training_scores)
] <- "high"

training_edges$fc_score <- binary_scores

# -----------------------------
# PREPARE TEST DATA
# -----------------------------

test_edges <- edge_features[
  is.na(edge_features$fc_score),
]

test_edges <- test_edges[
  ,
  -match(
    c(
      "fc_score",
      "fc_score_binned"
    ),
    colnames(test_edges)
  )
]

# -----------------------------
# FEATURE RENAMING
# -----------------------------

feature_names <- c(
  "cfreq",
  "k562i",
  "regde",
  "tarde",
  "regsur",
  "tarsurv",
  "corLIHC",
  "corGTEx",
  "corpLIHC",
  "corpGTEx",
  "regtlrbp",
  "tartlrbp",
  "regcan",
  "tarcan",
  "reggob",
  "targob",
  "regclit",
  "tarclit",
  "regcdeg",
  "tarcdeg",
  "eclip_reg",
  "eclip_targ",
  "fc_score"
)

colnames(training_edges) <- feature_names
colnames(test_edges) <- feature_names[-length(feature_names)]

# -----------------------------
# TRAIN CLASSIFICATION MODEL
# -----------------------------

training_edges$k562i <- as.numeric(
  training_edges$k562i
)

classification_model <- rpart(
  fc_score ~ .,
  training_edges
)

# Model visualization.
rpart.plot(
  classification_model,
  type = 5,
  tweak = 2,
  Margin = 0
)

# -----------------------------
# MODEL EVALUATION
# -----------------------------

cross_validation_predictions <- predict(
  classification_model,
  type = "class"
)

true_labels <- as.factor(
  training_edges$fc_score
)

confusionMatrix(
  cross_validation_predictions,
  true_labels
)

# -----------------------------
# PRUNE MODEL
# -----------------------------

pruned_model <- prune(
  classification_model,
  cp = 0.01414
)

rpart.plot(
  pruned_model,
  type = 5,
  tweak = 1.1,
  cex = 0.7
)

# -----------------------------
# PREDICT UNLABELED INTERACTIONS
# -----------------------------

test_edges$k562i <- as.numeric(
  test_edges$k562i
)

predicted_scores <- predict(
  classification_model,
  test_edges,
  type = "class"
) %>%
  as.data.frame()

test_edges$fc_score <- unname(
  as.matrix(predicted_scores)
)

# -----------------------------
# MERGE TRAINING + TEST DATA
# -----------------------------

merged_edges <- rbind(
  training_edges,
  test_edges
)

original_training <- edge_features[
  !is.na(edge_features$fc_score),
  -c(1, 2)
]

original_test <- edge_features[
  is.na(edge_features$fc_score),
  -c(1, 2)
]

merged_edges$fc_score_binned <-
  merged_edges$fc_score

merged_edges$fc_score <- c(
  original_training$fc_score,
  original_test$fc_score
)

saveRDS(
  merged_edges,
  OUTPUT_FINAL_EDGES
)

# -----------------------------
# BUILD MULTIFEATURE SCORE TABLE
# -----------------------------

scored_edges <- readRDS(
  OUTPUT_FINAL_EDGES
)

scored_edges$k562i[
  scored_edges$k562i == 0
] <- NaN

validation_table <- data.frame(
  GRN_frequency =
    abs(scored_edges$cfreq) /
    max(abs(scored_edges$cfreq)),
  
  in_k562 =
    scored_edges$k562i /
    max(scored_edges$k562i, na.rm = TRUE),
  
  DEG_LIHC =
    rowMeans(
      data.frame(
        scored_edges$regde,
        scored_edges$tarde
      ),
      na.rm = TRUE
    ),
  
  alter_survival =
    rowMeans(
      data.frame(
        scored_edges$regsur,
        scored_edges$tarsurv
      ),
      na.rm = TRUE
    ),
  
  coexpression_change =
    Mclust(
      abs(
        scored_edges$corLIHC -
          scored_edges$corGTEx
      ),
      N_CLUSTERS
    )$classification / N_CLUSTERS,
  
  regulator_eCLIP =
    cluster_with_missing(
      scored_edges$eclip_reg,
      2
    ) / 2,
  
  target_eCLIP =
    cluster_with_missing(
      scored_edges$eclip_targ,
      2
    ) / 2,
  
  fc_score =
    scored_edges$fc_score_binned
)

rownames(validation_table) <-
  rownames(scored_edges)

# -----------------------------
# COMPUTE FINAL RANKING SCORE
# -----------------------------

ranking_scores <- apply(
  validation_table[
    ,
    -length(validation_table)
  ],
  1,
  function(x) mean(
    na.omit(x)
  )
) %>%
  as.data.frame()

ranking_scores$score <-
  scored_edges$fc_score_binned %>%
  as.character()

colnames(ranking_scores) <- c(
  "rank_score",
  "fc_score"
)

# -----------------------------
# RANK INTERACTIONS
# -----------------------------

ranked_validation_table <-
  validation_table[
    order(
      ranking_scores$rank_score,
      decreasing = TRUE
    ),
  ][1:TOP_INTERACTIONS, ]

ranked_network <- validation_table[
  order(
    ranking_scores$rank_score,
    decreasing = TRUE
  ),
]

saveRDS(
  ranked_network,
  OUTPUT_RANKED_NETWORK
)

# -----------------------------
# NETWORK VISUALIZATION
# -----------------------------

interaction_network <- do.call(
  rbind,
  strsplit(
    rownames(ranked_validation_table),
    "→"
  )
) %>%
  as.data.frame()

colnames(interaction_network) <- c(
  "from",
  "to"
)

network_plot <- ggplot(
  ggnetwork(
    interaction_network,
    arrow.gap = 0.02
  ),
  aes(
    x,
    y,
    xend = xend,
    yend = yend
  )
) +
  geom_edges(
    arrow = arrow(
      length = unit(4, "pt"),
      type = "closed"
    ),
    curvature = 0.3,
    size = 1
  ) +
  geom_nodes(
    size = 5,
    color = "#636363"
  ) +
  geom_nodetext(
    aes(label = vertex.names),
    fontface = "bold",
    size = 4
  ) +
  theme_blank()

ggsave(
  OUTPUT_NETWORK_VISUALIZATION,
  plot = network_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# -----------------------------
# SAVE TOP NETWORK
# -----------------------------

write.csv2(
  interaction_network[, 1:2],
  OUTPUT_TOP_NETWORK,
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# ENRICHMENT ANALYSIS
# -----------------------------

top_genes <- unique(
  unlist(
    strsplit(
      rownames(ranked_validation_table),
      "→"
    )
  )
)

rbp_annotation <- read_excel(
  RBP_ANNOTATION_FILE
)

background_genes <- rbp_annotation$gene_name

gene_entrez <- mapIds(
  org.Hs.eg.db,
  top_genes,
  "ENTREZID",
  "SYMBOL"
)

background_entrez <- mapIds(
  org.Hs.eg.db,
  background_genes,
  "ENTREZID",
  "SYMBOL"
)

disease_enrichment <- enrichDGN(
  gene_entrez %>% unname(),
  universe = background_entrez %>% unname()
)

barplot(
  disease_enrichment,
  showCategory = 20
)

# ============================================================
# END OF SCRIPT
# ============================================================