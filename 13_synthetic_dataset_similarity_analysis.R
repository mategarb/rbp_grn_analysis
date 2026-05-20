# ============================================================
# Script: 13_syncode_synthetic_dataset_similarity_analysis.R
# Purpose: Evaluate similarity between real perturbation-based
#          transcriptomic datasets and SYNCODE-generated
#          synthetic datasets.
#
# This script:
#   - preprocesses HepG2/K562 perturbation matrices,
#   - compares synthetic and real distributions,
#   - evaluates similarity across SYNCODE parameter settings,
#   - ranks synthetic datasets using summary statistics,
#   - visualizes synthetic-vs-real similarity,
#   - identifies optimal synthetic configurations.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files
HEPG2_FILE <- file.path(
  INPUT_DIR,
  "ymatrix_hepg2.csv"
)

K562_FILE <- file.path(
  INPUT_DIR,
  "ymatrix_k562.csv"
)

SYNCODE_DIR <- file.path(
  INPUT_DIR,
  "SYNCODE_results"
)

# Select dataset
CELL_LINE <- "K562"   # Options: "HepG2", "K562"

# Genes removed before comparison
GENES_TO_REMOVE <- c(
  "CELF1",
  "DDX24",
  "DDX28",
  "DNAJC21",
  "EWSR1",
  "HNRNPA0",
  "HNRNPC",
  "IGF2BP3",
  "PABPC1",
  "PARN",
  "RCC2",
  "RECQL",
  "RPS3A",
  "SF1",
  "TARDBP",
  "TFIP11"
)

# Output files
OUTPUT_HEATMAP_REAL <- file.path(
  OUTPUT_DIR,
  "real_expression_heatmap.svg"
)

OUTPUT_HEATMAP_SYNTHETIC <- file.path(
  OUTPUT_DIR,
  "best_synthetic_expression_heatmap.svg"
)

OUTPUT_SIMILARITY_BOXPLOT <- file.path(
  OUTPUT_DIR,
  "syncode_similarity_scores.svg"
)

OUTPUT_RANKING_PLOT <- file.path(
  OUTPUT_DIR,
  "syncode_ranking.svg"
)

# -----------------------------
# LIBRARIES
# -----------------------------
library(tidyverse)
library(pheatmap)
library(R.matlab)
library(colordistance)
library(stringr)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Scale values to 0-1 range.
scale_01 <- function(x) {
  
  (x - min(x)) /
    (max(x) - min(x))
}

# Compute summary-statistic distance.
compute_summary_distance <- function(
    real_data,
    synthetic_data
) {
  
  real_vector <- as.vector(
    as.matrix(real_data)
  )
  
  synthetic_vector <- as.vector(
    as.matrix(synthetic_data)
  )
  
  abs(
    summary(real_vector) -
      summary(synthetic_vector)
  )
}

# -----------------------------
# LOAD EXPRESSION DATA
# -----------------------------

# HepG2
hepg2_raw <- read.csv2(
  file = HEPG2_FILE,
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

# K562
k562_raw <- read.csv2(
  file = K562_FILE,
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

# Select dataset.
input_data <- switch(
  CELL_LINE,
  "HepG2" = hepg2_expression,
  "K562" = k562_expression
)

# -----------------------------
# REMOVE SELECTED GENES
# -----------------------------

rows_to_remove <- match(
  GENES_TO_REMOVE,
  rownames(input_data)
)

cols_to_remove <- sapply(
  GENES_TO_REMOVE,
  function(gene) {
    
    grep(
      paste0("Gene", gene),
      paste0("Gene", colnames(input_data)),
      fixed = TRUE
    )
  }
) %>%
  as.vector()

input_data <- input_data[
  -rows_to_remove,
  -cols_to_remove
]

# -----------------------------
# VISUALIZE REAL DATA
# -----------------------------

svg(
  OUTPUT_HEATMAP_REAL,
  width = 10,
  height = 8
)

pheatmap(
  scale(input_data),
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

dev.off()

# -----------------------------
# LOAD SYNCODE DATASETS
# -----------------------------

syncode_files <- list.files(
  SYNCODE_DIR
)

syncode_files <- syncode_files[
  str_detect(syncode_files, "Y_")
]

# -----------------------------
# COMPARE SYNTHETIC DATASETS
# -----------------------------

distance_scores <- numeric(
  length(syncode_files)
)

summary_distance_list <- list()

for (i in seq_along(syncode_files)) {
  
  synthetic_data <- readMat(
    file.path(
      SYNCODE_DIR,
      syncode_files[i]
    )
  )
  
  synthetic_expression <- synthetic_data$Yo %>%
    as.data.frame()
  
  # Standardize matrices.
  real_scaled <- scale(input_data)
  synthetic_scaled <- scale(synthetic_expression)
  
  # Summary-statistic distance.
  summary_distance_list[[i]] <-
    compute_summary_distance(
      real_scaled,
      synthetic_scaled
    )
  
  distance_scores[i] <-
    sum(summary_distance_list[[i]])
}

# -----------------------------
# IDENTIFY BEST SYNTHETIC DATASET
# -----------------------------

best_dataset_index <- which.min(
  distance_scores
)

best_dataset_file <- syncode_files[
  best_dataset_index
]

best_syncode <- readMat(
  file.path(
    SYNCODE_DIR,
    best_dataset_file
  )
)

best_syncode_expression <- best_syncode$Yo %>%
  as.data.frame()

# -----------------------------
# VISUALIZE BEST SYNTHETIC DATASET
# -----------------------------

svg(
  OUTPUT_HEATMAP_SYNTHETIC,
  width = 10,
  height = 8
)

pheatmap(
  scale(best_syncode_expression) %>%
    scale_01(),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = best_dataset_file
)

dev.off()

# -----------------------------
# EXTRACT SYNCODE PARAMETERS
# -----------------------------

syncode_parameters <- str_match(
  syncode_files,
  "Y_SYNCODE_\\s*(.*?)\\_ITER"
)[, 2]

# -----------------------------
# SIMILARITY SCORE BOXPLOT
# -----------------------------

similarity_data <- data.frame(
  parameter = syncode_parameters,
  score = distance_scores
)

p_similarity <- ggplot(
  similarity_data,
  aes(
    x = parameter,
    y = score
  )
) +
  geom_boxplot(
    outlier.colour = "red",
    outlier.shape = 8,
    outlier.size = 4
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red"
  ) +
  theme_minimal() +
  ylab("Similarity score") +
  xlab("SYNCODE parameter")

ggsave(
  OUTPUT_SIMILARITY_BOXPLOT,
  plot = p_similarity,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# -----------------------------
# RANK SYNTHETIC DATASETS
# -----------------------------

ranking_group <- rep(
  "remaining datasets",
  length(syncode_files)
)

top_indices <- order(distance_scores)[1:5]

ranking_group[top_indices] <-
  gsub(
    "_S=3.mat",
    "",
    gsub(
      "Y_SYNCODE_",
      "",
      syncode_files
    )
  )[top_indices]

ranking_df <- data.frame(
  syncode = gsub(
    "_S=3.mat",
    "",
    gsub(
      "Y_SYNCODE_",
      "",
      syncode_files
    )
  ),
  score = distance_scores,
  group = ranking_group
)

ranking_df <- ranking_df[
  order(ranking_df$score),
]

# -----------------------------
# VISUALIZE RANKING
# -----------------------------

p_ranking <- ggplot(
  ranking_df,
  aes(
    x = reorder(syncode, score),
    y = score,
    fill = reorder(group, score)
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "darkred",
      "brown1",
      "darkorange",
      "darkgoldenrod1",
      "gold",
      "slategray"
    ),
    name = ""
  ) +
  ylab("Summary statistics score") +
  xlab("SYNCODE dataset") +
  ggtitle(CELL_LINE) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  OUTPUT_RANKING_PLOT,
  plot = p_ranking,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ============================================================
# END OF SCRIPT
# ============================================================