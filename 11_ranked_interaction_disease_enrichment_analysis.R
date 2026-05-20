# ============================================================
# Script: 11_ranked_interaction_disease_enrichment_analysis.R
# Purpose: Identify disease-enriched regions within ranked
#          interaction networks using incremental thresholds
#          and sliding-window enrichment analyses.
#
# This script:
#   - evaluates disease enrichment across ranked interactions,
#   - prioritizes liver/hepatocellular disease-associated links,
#   - applies harmonic mean p-value aggregation,
#   - identifies significant interaction regions,
#   - visualizes enrichment trajectories,
#   - and exports prioritized genes.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files
RANKED_SCORES_FILE <- file.path(
  INPUT_DIR,
  "net_scores_val.rds"
)

RBP_ANNOTATION_FILE <- file.path(
  INPUT_DIR,
  "rbps_info",
  "210329_Table_S1_hRBP_list.xlsx"
)

# Output files
OUTPUT_ENRICHMENT_PLOT <- file.path(
  OUTPUT_DIR,
  "ranked_interaction_disease_enrichment.svg"
)

OUTPUT_TOP_GENES <- file.path(
  OUTPUT_DIR,
  "top_genes.txt"
)

OUTPUT_ENRICHMENT_BARPLOT <- file.path(
  OUTPUT_DIR,
  "final_disease_enrichment.svg"
)

# Parameters
SLIDING_WINDOW_SIZE <- 30
SIGNIFICANCE_THRESHOLD <- 0.05

# -----------------------------
# LIBRARIES
# -----------------------------
library(gprofiler2)
library(readxl)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(harmonicmeanp)
library(ggpubr)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Convert gene symbols to Entrez IDs.
convert_to_entrez <- function(genes) {
  
  mapIds(
    org.Hs.eg.db,
    genes,
    "ENTREZID",
    "SYMBOL"
  )
}

# Extract liver/hepatocellular enrichment terms.
extract_liver_terms <- function(enrichment_results) {
  
  if (nrow(enrichment_results) == 0) {
    return(NULL)
  }
  
  liver_indices <- unique(c(
    grep("liver", enrichment_results$Description),
    grep("Liver", enrichment_results$Description),
    grep("hepato", enrichment_results$Description),
    grep("Hepato", enrichment_results$Description)
  ))
  
  enrichment_results[liver_indices, ]
}

# Run disease enrichment analysis.
run_disease_enrichment <- function(
    genes,
    background_genes
) {
  
  gene_entrez <- convert_to_entrez(genes)
  
  enrichDGN(
    unname(gene_entrez),
    universe = unname(background_genes)
  )
}

# Aggregate p-values using harmonic mean p-value.
harmonic_mean_pvalues <- function(pvalue_list) {
  
  aggregated_pvalues <- numeric(length(pvalue_list))
  
  for (i in seq_along(pvalue_list)) {
    
    if (
      length(pvalue_list[[i]]) <= 1 ||
      pvalue_list[[i]][1] == 0
    ) {
      
      aggregated_pvalues[i] <- 1
      
    } else {
      
      aggregated_pvalues[i] <- p.hmp(
        pvalue_list[[i]],
        L = length(pvalue_list[[i]])
      ) %>%
        unname()
    }
  }
  
  aggregated_pvalues
}

# Extract genes from ranked interaction rows.
extract_genes_from_interactions <- function(interaction_rows) {
  
  unique(
    unlist(
      strsplit(
        interaction_rows,
        "→"
      )
    )
  )
}

# -----------------------------
# LOAD DATA
# -----------------------------

ranked_scores <- readRDS(
  RANKED_SCORES_FILE
)

rbp_annotation <- read_excel(
  RBP_ANNOTATION_FILE
)

background_genes <- rbp_annotation$gene_name

background_entrez <- convert_to_entrez(
  background_genes
)

thresholds <- seq_len(
  nrow(ranked_scores)
)

# -----------------------------
# INCREMENTAL THRESHOLD ANALYSIS
# -----------------------------

incremental_pvalues <- list()
incremental_adjusted <- list()

for (i in thresholds) {
  
  selected_genes <- extract_genes_from_interactions(
    rownames(
      ranked_scores[1:i, ]
    )
  )
  
  enrichment <- run_disease_enrichment(
    selected_genes,
    background_entrez
  )
  
  enrichment_results <- enrichment@result
  
  liver_terms <- extract_liver_terms(
    enrichment_results
  )
  
  if (
    is.null(liver_terms) ||
    nrow(liver_terms) == 0
  ) {
    
    incremental_pvalues[[i]] <- 0
    incremental_adjusted[[i]] <- 0
    
  } else {
    
    incremental_pvalues[[i]] <-
      liver_terms$pvalue
    
    names(incremental_pvalues[[i]]) <-
      liver_terms$Description
    
    incremental_adjusted[[i]] <-
      liver_terms$p.adjust
    
    names(incremental_adjusted[[i]]) <-
      liver_terms$Description
  }
  
  print(i)
}

# -----------------------------
# SLIDING-WINDOW ANALYSIS
# -----------------------------

window_pvalues <- list()
window_adjusted <- list()

for (
  i in seq_len(
    length(thresholds) -
    SLIDING_WINDOW_SIZE
  )
) {
  
  selected_genes <- extract_genes_from_interactions(
    rownames(
      ranked_scores[
        i:(i + SLIDING_WINDOW_SIZE),
      ]
    )
  )
  
  enrichment <- run_disease_enrichment(
    selected_genes,
    background_entrez
  )
  
  enrichment_results <- enrichment@result
  
  liver_terms <- extract_liver_terms(
    enrichment_results
  )
  
  if (
    is.null(liver_terms) ||
    nrow(liver_terms) == 0
  ) {
    
    window_pvalues[[i]] <- 0
    window_adjusted[[i]] <- 0
    
  } else {
    
    window_pvalues[[i]] <-
      liver_terms$pvalue
    
    names(window_pvalues[[i]]) <-
      liver_terms$Description
    
    window_adjusted[[i]] <-
      liver_terms$p.adjust
    
    names(window_adjusted[[i]]) <-
      liver_terms$Description
  }
  
  print(i)
}

# -----------------------------
# HARMONIC-MEAN P-VALUE SUMMARY
# -----------------------------

incremental_summary <- harmonic_mean_pvalues(
  incremental_adjusted
)

window_summary <- harmonic_mean_pvalues(
  window_adjusted
)

# -----------------------------
# IDENTIFY SIGNIFICANT REGIONS
# -----------------------------

best_incremental_threshold <- which.min(
  incremental_summary
)

top_interactions <- seq_len(
  best_incremental_threshold
)

significant_windows <- which(
  window_summary <= SIGNIFICANCE_THRESHOLD
)

window_breaks <- which(
  diff(significant_windows) != 1
)

window_start <- numeric()
window_end <- numeric()

window_start[1] <- 1

for (i in seq_along(window_breaks)) {
  
  window_end[i] <-
    significant_windows[window_breaks[i]]
  
  window_start[i + 1] <-
    significant_windows[
      window_breaks[i] + 1
    ]
}

window_end[length(window_breaks) + 1] <-
  significant_windows[
    length(significant_windows)
  ]

significant_regions <- data.frame(
  start = window_start,
  end = window_end
)

significant_regions_extended <-
  significant_regions

significant_regions_extended$end <-
  significant_regions_extended$end +
  SLIDING_WINDOW_SIZE

# -----------------------------
# VISUALIZATION
# -----------------------------

incremental_plot_data <- data.frame(
  threshold = thresholds[
    seq_along(incremental_summary)
  ],
  pvalue = incremental_summary
)

p_incremental <- ggplot(
  incremental_plot_data,
  aes(
    x = threshold,
    y = pvalue
  )
) +
  geom_bar(
    stat = "identity",
    fill = "#7EAFC7"
  ) +
  geom_hline(
    yintercept = SIGNIFICANCE_THRESHOLD,
    linetype = "dashed",
    color = "#FE3200"
  ) +
  geom_vline(
    xintercept = best_incremental_threshold,
    linetype = "dashed",
    color = "#FE3200"
  ) +
  theme_minimal() +
  xlab("Number of top interactions") +
  ylab("P value")

window_plot_data <- data.frame(
  threshold = thresholds[
    seq_along(window_summary)
  ],
  pvalue = window_summary
)

p_window <- ggplot(
  window_plot_data,
  aes(
    x = threshold,
    y = pvalue
  )
) +
  geom_bar(
    stat = "identity",
    fill = "#7EAFC7"
  ) +
  geom_hline(
    yintercept = SIGNIFICANCE_THRESHOLD,
    linetype = "dashed",
    color = "#FE3200"
  ) +
  geom_rect(
    data = significant_regions,
    inherit.aes = FALSE,
    aes(
      xmin = start,
      xmax = end,
      ymin = min(window_plot_data$pvalue),
      ymax = max(window_plot_data$pvalue)
    ),
    fill = "#FE3200",
    alpha = 0.4,
    color = NA
  ) +
  theme_minimal() +
  xlab(
    paste0(
      "Sliding window size = ",
      SLIDING_WINDOW_SIZE
    )
  ) +
  ylab("P value")

combined_plot <- ggarrange(
  p_incremental,
  p_window,
  nrow = 2
)

ggsave(
  OUTPUT_ENRICHMENT_PLOT,
  combined_plot,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

# -----------------------------
# FINAL ENRICHMENT ANALYSIS
# -----------------------------

window_interactions <- c()

for (
  i in seq_len(
    nrow(significant_regions_extended)
  )
) {
  
  window_interactions <- c(
    window_interactions,
    significant_regions_extended$start[i]:
      significant_regions_extended$end[i]
  )
}

selected_indices <- union(
  top_interactions,
  unique(window_interactions)
)

selected_genes <- extract_genes_from_interactions(
  rownames(
    ranked_scores[selected_indices, ]
  )
)

final_enrichment <- run_disease_enrichment(
  selected_genes,
  background_entrez
)

# -----------------------------
# SAVE FINAL OUTPUTS
# -----------------------------

svg(
  OUTPUT_ENRICHMENT_BARPLOT,
  width = 10,
  height = 8
)

barplot(
  final_enrichment,
  showCategory = 20
)

dev.off()

write.table(
  paste(
    sort(unique(selected_genes)),
    collapse = ", "
  ),
  OUTPUT_TOP_GENES,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# ============================================================
# END OF SCRIPT
# ============================================================