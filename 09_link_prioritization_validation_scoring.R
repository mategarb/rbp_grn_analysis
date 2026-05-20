# ============================================================
# Script: 09_link_prioritization_validation_scoring.R
# Purpose: Prioritize candidate RBP-RBP interactions using multi-source
#          validation evidence, including LIHC/GTEx differential expression,
#          survival association, co-expression change, RBP/cancer annotations,
#          FunCoup scores, eCLIP support, and disease enrichment.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files
VALIDATION_DIR <- file.path(INPUT_DIR, "validation_UCSC_LIHC_GTEx/")
RBP_INFO_DIR <- file.path(INPUT_DIR, "rbps_info/")
ECLIP_DIR <- file.path(INPUT_DIR, "eCLIP_data/")
RESULTS_DIR <- file.path(INPUT_DIR, "results/")

GTEX_PHENOTYPE_FILE <- file.path(VALIDATION_DIR, "GTEX_phenotype.gz")
LIHC_SURVIVAL_FILE <- file.path(VALIDATION_DIR, "survival_LIHC_survival.txt")
GENE_PROBEMAP_FILE <- file.path(VALIDATION_DIR, "probeMap_gencode.v23.annotation.gene.probemap")
LIVER_EXPRESSION_FILE <- file.path(RESULTS_DIR, "liver_tcga_gtex.rds")
ALL_INTERACTIONS_FILE <- file.path(RESULTS_DIR, "all_inters_hepg2.rds")

RBP_LIST_FILE <- file.path(RBP_INFO_DIR, "210329_Table_S1_hRBP_list.xlsx")
CANCER_LITERATURE_FILE <- file.path(RBP_INFO_DIR, "TheNumberOfCancerRelevantLiteraturesOfAllRBPs.csv")
CANCER_DEG_FILE <- file.path(RBP_INFO_DIR, "TheNumberOfDifferentiallyExpressedRBPsOfAllCancerTypes.csv")
FUNCOUP_NETWORK_FILE <- file.path(INPUT_DIR, "FC5.0_H.sapiens_full.gz")
GENCODE_ANNOTATION_FILE <- file.path(INPUT_DIR, "gencode.annotation.gtf.gz")
ECLIP_METADATA_FILE <- file.path(ECLIP_DIR, "metadata.tsv")

# Output files
OUTPUT_EDGE_SUMMARY_RDS <- file.path(OUTPUT_DIR, "edges_interaction_validation_summary.rds")
OUTPUT_FINAL_EDGE_TABLE_RDS <- file.path(OUTPUT_DIR, "edges_final_hepg2.rds")
OUTPUT_ECLIP_VALIDATED_RDS <- file.path(OUTPUT_DIR, "interactions_hepg2_eclip.rds")
OUTPUT_RANKED_LINKS_CSV <- file.path(OUTPUT_DIR, "all_ranked_cGRN_links.csv")
OUTPUT_TOP_INTERACTIONS_CSV <- file.path(OUTPUT_DIR, "top_interactions.csv")
OUTPUT_VALIDATION_HEATMAP <- file.path(OUTPUT_DIR, "validation_feature_table.svg")
OUTPUT_FUNCOUP_NETWORK_FIGURE <- file.path(OUTPUT_DIR, "network_funcoup_score.svg")
OUTPUT_DISEASE_ENRICHMENT_FIGURE <- file.path(OUTPUT_DIR, "disease_enrichment_ranking.svg")

# Analysis parameters
TOP_N_INTERACTIONS <- 27
TOP_N_FOR_WORDCLOUD <- 130
N_MCLUST_GROUPS <- 3
ECLIP_FDR_THRESHOLD <- 0.05
DISEASE_WINDOW_SIZE <- 30
DISEASE_PVALUE_THRESHOLD <- 0.05

# -----------------------------
# LIBRARIES
# -----------------------------
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(scales)
library(classInt)
library(rpart)
library(rpart.plot)
library(readxl)
library(rtracklayer)
library(GenomicRanges)
library(RCAS)
library(org.Hs.eg.db)
library(caret)
library(mclust)
library(DOSE)
library(enrichplot)
library(ggnetwork)
library(igraph)
library(clusterProfiler)
library(harmonicmeanp)
library(ggpmisc)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Discretize p-values into significance-score categories.
discretize_pvalue <- function(p) {
  case_when(
    is.nan(p) ~ 0,
    p > 0.1 ~ 0,
    p <= 0.1 & p > 0.05 ~ 1,
    p <= 0.05 & p > 0.01 ~ 2,
    p <= 0.01 & p > 0.001 ~ 3,
    p <= 0.001 ~ 4,
    TRUE ~ NA_real_
  )
}

# Convert discretized p-value scores into symbols.
pvalue_score_to_symbol <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x == 0 ~ "ns",
    x == 1 ~ "*",
    x == 2 ~ "**",
    x == 3 ~ "***",
    x == 4 ~ "****",
    TRUE ~ "NA"
  )
}

# Discretize absolute correlations.
discretize_correlation <- function(x) {
  x <- abs(round(x, digits = 1))
  case_when(
    is.nan(x) ~ 0,
    x >= 0.7 ~ 4,
    x < 0.7 & x >= 0.5 ~ 3,
    x < 0.5 & x >= 0.3 ~ 2,
    x < 0.3 & x > 0 ~ 1,
    TRUE ~ 0
  )
}

# Harmonize symbols used by different resources.
to_expression_symbol <- function(x) {
  x[x == "ATP5F1C"] <- "ATP5C1"
  x[x == "RACK1"] <- "GNB2L1"
  x
}

to_public_symbol <- function(x) {
  x[x == "ATP5C1"] <- "ATP5F1C"
  x[x == "GNB2L1"] <- "RACK1"
  x
}

# Remove repeated expression values that dominate a gene distribution.
filter_repeated_expression_values <- function(df, expression_columns) {
  keep <- rep(TRUE, nrow(df))
  
  for (col in expression_columns) {
    value_table <- table(df[[col]])
    if (length(which(value_table > 10)) >= 2) {
      repeated_values <- names((sort(value_table) / sum(value_table))[sort(value_table) / sum(value_table) >= 0.05])
      keep <- keep & !(as.numeric(df[[col]]) %in% as.numeric(repeated_values))
    }
  }
  
  df[keep, , drop = FALSE]
}

# Run DEG, survival, and co-expression statistics for one interaction.
analyze_interaction_expression <- function(gene1, gene2, expression_matrix, survival_table) {
  genes <- to_expression_symbol(c(gene1, gene2))
  
  gene_stats <- lapply(genes, function(gene) {
    gene_col <- match(gene, colnames(expression_matrix))
    if (is.na(gene_col)) {
      return(list(deg_p = NA_real_, survival_p = NA_real_))
    }
    
    gene_expr <- as.data.frame(expression_matrix[gene_col])
    gene_expr$id <- substr(rownames(gene_expr), 1, 4)
    gene_expr$id[gene_expr$id == "TCGA"] <- "LIHC"
    gene_expr$id[gene_expr$id == "GTEX"] <- "control"
    gene_expr <- filter_repeated_expression_values(gene_expr, names(gene_expr)[1])
    
    deg_test <- t.test(gene_expr[gene_expr$id == "LIHC", 1], gene_expr[gene_expr$id == "control", 1])
    
    survival_gene <- survival_table[match(rownames(gene_expr), survival_table$sample) %>% na.omit()] %>% as.data.frame()
    survival_expr <- gene_expr[match(survival_table$sample, rownames(gene_expr)) %>% na.omit(), 1]
    survival_gene$ge_status <- ifelse(survival_expr <= median(survival_expr), "low", "high")
    
    fit <- survfit(Surv(OS.time, OS) ~ ge_status, data = survival_gene)
    survival_p <- surv_pvalue(fit, survival_gene)$pval
    
    list(deg_p = deg_test$p.value, survival_p = survival_p)
  })
  
  pair_cols <- match(genes, colnames(expression_matrix))
  if (any(is.na(pair_cols))) {
    return(c(gene_stats[[1]]$deg_p, gene_stats[[2]]$deg_p, gene_stats[[1]]$survival_p, gene_stats[[2]]$survival_p, NA, NA, NA, NA))
  }
  
  pair_expr <- expression_matrix[pair_cols]
  colnames(pair_expr) <- c("gene1", "gene2")
  pair_expr <- filter_repeated_expression_values(pair_expr, c("gene1", "gene2"))
  pair_expr$id <- substr(rownames(pair_expr), 1, 4)
  
  tcga_cor <- cor.test(pair_expr$gene1[pair_expr$id == "TCGA"], pair_expr$gene2[pair_expr$id == "TCGA"])
  
  if (length(pair_expr$gene1[pair_expr$id == "GTEX"]) == 0) {
    gtex_cor_p <- NaN
    gtex_cor <- NaN
  } else {
    gtex_cor_test <- cor.test(pair_expr$gene1[pair_expr$id == "GTEX"], pair_expr$gene2[pair_expr$id == "GTEX"])
    gtex_cor_p <- gtex_cor_test$p.value
    gtex_cor <- as.numeric(gtex_cor_test$estimate)
  }
  
  c(
    gene_stats[[1]]$deg_p,
    gene_stats[[2]]$deg_p,
    gene_stats[[1]]$survival_p,
    gene_stats[[2]]$survival_p,
    as.numeric(tcga_cor$estimate),
    gtex_cor,
    tcga_cor$p.value,
    gtex_cor_p
  )
}

# Add RBP and cancer annotations to interaction table.
add_rbp_cancer_annotations <- function(edges, rbp_table, cancer_lit, cancer_deg) {
  edges$regulator_times_listed_as_rbp <- rbp_table$Times_Listed_as_RBP[match(edges$from, rbp_table$gene_name)]
  edges$target_times_listed_as_rbp <- rbp_table$Times_Listed_as_RBP[match(edges$to, rbp_table$gene_name)]
  edges$regulator_canonical <- rbp_table$`canonical/non_canonical`[match(edges$from, rbp_table$gene_name)]
  edges$target_canonical <- rbp_table$`canonical/non_canonical`[match(edges$to, rbp_table$gene_name)]
  edges$regulator_go_binding <- rbp_table$Gene_Ontology_RNA_Binding[match(edges$from, rbp_table$gene_name)]
  edges$target_go_binding <- rbp_table$Gene_Ontology_RNA_Binding[match(edges$to, rbp_table$gene_name)]
  edges$regulator_cancer_lit <- cancer_lit$Number.of.Literatures[match(edges$from, cancer_lit$Gene.Symbol)]
  edges$target_cancer_lit <- cancer_lit$Number.of.Literatures[match(edges$to, cancer_lit$Gene.Symbol)]
  edges$regulator_cancer_deg <- cancer_deg$Number.of.cancers[match(edges$from, cancer_deg$Gene.Symbol)]
  edges$target_cancer_deg <- cancer_deg$Number.of.cancers[match(edges$to, cancer_deg$Gene.Symbol)]
  
  edges
}

# Attach FunCoup score labels using known FunCoup interactions as training labels.
add_funcoup_scores <- function(edges, funcoup_network) {
  genes_ensembl <- sapply(seq_len(nrow(edges)), function(i) {
    symbols <- mapIds(org.Hs.eg.db, keys = c(edges$from[i], edges$to[i]), keytype = "SYMBOL", column = "ENSEMBL")
    paste(unname(symbols), collapse = "-")
  })
  
  funcoup_forward <- paste0(funcoup_network$`2:Gene1`, "-", funcoup_network$`3:Gene2`)
  funcoup_reverse <- paste0(funcoup_network$`3:Gene2`, "-", funcoup_network$`2:Gene1`)
  
  matched_indices <- sort(c(which(!is.na(match(funcoup_forward, genes_ensembl))), which(!is.na(match(funcoup_reverse, genes_ensembl)))))
  funcoup_subset <- funcoup_network[matched_indices, ]
  
  funcoup_subset$`2:Gene1` <- mapIds(org.Hs.eg.db, keys = funcoup_subset$`2:Gene1`, keytype = "ENSEMBL", column = "SYMBOL") %>% unname()
  funcoup_subset$`3:Gene2` <- mapIds(org.Hs.eg.db, keys = funcoup_subset$`3:Gene2`, keytype = "ENSEMBL", column = "SYMBOL") %>% unname()
  
  funcoup_match <- rep(NA_integer_, nrow(edges))
  for (i in seq_len(nrow(edges))) {
    query_pair <- paste0(edges$from[i], "-", edges$to[i])
    for (j in seq_len(nrow(funcoup_subset))) {
      fc_pair <- paste0(funcoup_subset$`2:Gene1`[j], "-", funcoup_subset$`3:Gene2`[j])
      fc_pair_rev <- paste0(funcoup_subset$`3:Gene2`[j], "-", funcoup_subset$`2:Gene1`[j])
      if (query_pair == fc_pair || query_pair == fc_pair_rev) {
        funcoup_match[i] <- j
      }
    }
  }
  
  edges$fc_score <- funcoup_subset$`#0:PFC`[funcoup_match]
  edges
}

# Train decision tree to impute binned FunCoup classes for unmatched links.
impute_funcoup_classes <- function(edges) {
  model_edges <- edges
  model_edges$weight <- abs(model_edges$weight)
  
  train_edges <- model_edges[!is.na(model_edges$fc_score), -c(1, 2)]
  test_edges <- model_edges[is.na(model_edges$fc_score), -c(1, 2)]
  rownames(train_edges) <- paste0(model_edges$from, "-", model_edges$to)[!is.na(model_edges$fc_score)]
  rownames(test_edges) <- paste0(model_edges$from, "-", model_edges$to)[is.na(model_edges$fc_score)]
  
  breaks <- classIntervals(train_edges$fc_score, 2, style = "equal")
  breaks$brks[1] <- 0
  train_edges$fc_score <- cut(breaks$var, breaks$brks, labels = c("low", "high"))
  
  new_names <- c(
    "cfreq", "k562i", "regde", "tarde", "regsur", "tarsurv", "corLIHC", "corGTEx",
    "corpLIHC", "corpGTEx", "regtlrbp", "tartlrbp", "regcan", "tarcan", "reggob", "targob",
    "regclit", "tarclit", "regcdeg", "tarcdeg", "fc_score"
  )
  colnames(train_edges) <- new_names
  colnames(test_edges) <- new_names
  
  train_edges$k562i <- as.numeric(as.matrix(train_edges$k562i))
  train_edges$cfreq <- abs(train_edges$cfreq)
  
  decision_tree <- rpart(fc_score ~ ., train_edges, method = "class", minsplit = 5, minbucket = 2)
  
  cv_prediction <- predict(decision_tree, type = "class")
  confusion <- confusionMatrix(cv_prediction, train_edges$fc_score %>% as.factor())
  
  test_edges$k562i <- as.numeric(as.matrix(test_edges$k562i))
  test_prediction <- predict(decision_tree, test_edges, type = "class") %>% as.data.frame()
  test_edges$fc_score <- unname(as.matrix(test_prediction))
  
  merged <- rbind(train_edges, test_edges)
  original_train <- model_edges[!is.na(model_edges$fc_score), -c(1, 2)]
  original_test <- model_edges[is.na(model_edges$fc_score), -c(1, 2)]
  merged$fc_score_binned <- merged$fc_score
  merged$fc_score <- c(original_train$fc_score, original_test$fc_score)
  
  list(edges = merged, model = decision_tree, confusion = confusion)
}

# Extract eCLIP support for a regulator-target pair.
get_eclip_pair_support <- function(regulator, target, metadata, annotation, eclip_dir, fdr_threshold = 0.05) {
  meta_sub <- metadata[gsub("-human", "", metadata$`Experiment target`) == regulator, ]
  
  if (nrow(meta_sub) == 0) return(NaN)
  
  file_ids <- meta_sub$`File accession`[meta_sub$`File format` == "bigBed narrowPeak"]
  if (length(file_ids) == 0) return(NaN)
  
  hit_count <- 0
  for (file_id in file_ids) {
    peak_file <- file.path(eclip_dir, paste0(file_id, ".bigBed"))
    if (!file.exists(peak_file)) next
    
    peak_data <- import.bb(peak_file)
    names(mcols(peak_data)) <- c("name", "score", "signalValue", "pValue", "qValue", "peak")
    
    overlaps <- as.data.table(queryGff(queryRegions = peak_data, gffData = annotation))
    overlaps <- overlaps[gene_type == "protein_coding" & type == "gene"]
    overlaps <- overlaps[p.adjust(10^(-as.numeric(query_pValue)), method = "fdr") <= fdr_threshold]
    
    hit_count <- hit_count + length(which(overlaps$gene_name == target))
  }
  
  hit_count
}

# Cluster numeric vector while preserving NA and zero values.
cluster_with_na_and_zero <- function(x, n_groups) {
  non_na <- which(!is.na(x))
  zero_idx <- which(x == 0)
  output <- rep(NA_real_, length(x))
  
  to_cluster <- setdiff(non_na, zero_idx)
  if (length(to_cluster) > 1) {
    output[to_cluster] <- Mclust(as.numeric(x[to_cluster]), n_groups)$classification
  }
  output[zero_idx] <- 0
  output
}

# Build normalized validation feature matrix and ranking score.
build_validation_score_matrix <- function(edges, n_mclust_groups = 3) {
  edges$k562i[edges$k562i == 0] <- NaN
  
  feature_df <- data.frame(
    GRN_frequency = abs(edges$cfreq) / max(abs(edges$cfreq), na.rm = TRUE),
    in_k562 = edges$k562i / max(edges$k562i, na.rm = TRUE),
    DEG_LIHC = rowMeans(data.frame(edges$regde, edges$tarde), na.rm = TRUE) / max(rowMeans(data.frame(edges$regde, edges$tarde), na.rm = TRUE), na.rm = TRUE),
    alter_survival = rowMeans(data.frame(edges$regsur, edges$tarsurv), na.rm = TRUE) / max(rowMeans(data.frame(edges$regsur, edges$tarsurv), na.rm = TRUE), na.rm = TRUE),
    coexpression_change = Mclust(abs(edges$corLIHC - edges$corGTEx), n_mclust_groups)$classification / n_mclust_groups,
    times_listed_as_RBP = Mclust(rowMeans(data.frame(edges$regtlrbp, edges$tartlrbp), na.rm = TRUE), n_mclust_groups)$classification / n_mclust_groups,
    in_literature_as_cancer = cluster_with_na_and_zero(rowMeans(data.frame(edges$regclit, edges$tarclit), na.rm = TRUE), 2) / 2,
    in_cancers_as_DEG = cluster_with_na_and_zero(rowMeans(data.frame(edges$regcdeg, edges$tarcdeg), na.rm = TRUE), 2) / 2,
    RBD = rowMeans(data.frame(edges$regcan, edges$tarcan), na.rm = TRUE) / max(rowMeans(data.frame(edges$regcan, edges$tarcan), na.rm = TRUE), na.rm = TRUE),
    regulator_eCLIP = cluster_with_na_and_zero(edges$eclip_reg, 2) / 2,
    target_eCLIP = cluster_with_na_and_zero(edges$eclip_targ, 2) / 2,
    fc_score = edges$fc_score_binned
  )
  
  rownames(feature_df) <- gsub("-", "→", rownames(edges))
  rank_score <- apply(feature_df[, -ncol(feature_df)], 1, function(x) mean(na.omit(x)))
  
  list(features = feature_df, rank_score = rank_score)
}

# Scale helper for ternary/summary plots.
scale_values <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Harmonic-mean p-value summary over list elements.
pval_mean <- function(pval_list) {
  sapply(pval_list, function(x) {
    if (length(x) <= 1 || x[1] == 0) {
      1
    } else {
      p.hmp(x, L = length(x)) %>% unname()
    }
  })
}

# Plot disease-enrichment signal across ranked interactions.
plot_rank_disease_enrichment <- function(ranked_features, rbp_background, window_size = 30) {
  thresholds <- seq_len(nrow(ranked_features))
  top_pvalues <- list()
  window_pvalues <- list()
  
  for (i in thresholds) {
    genes <- unique(unlist(strsplit(rownames(ranked_features[1:i, , drop = FALSE]), "→")))
    gene_entrez <- mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
    enrichment <- enrichDGN(unname(gene_entrez), universe = unname(rbp_background))
    result <- enrichment@result
    
    if (nrow(result) == 0) {
      top_pvalues[[i]] <- 0
    } else {
      liver_hits <- unique(c(grep("liver", result$Description), grep("Liver", result$Description), grep("hepato", result$Description), grep("Hepato", result$Description)))
      top_pvalues[[i]] <- result$p.adjust[liver_hits]
      names(top_pvalues[[i]]) <- result$Description[liver_hits]
    }
  }
  
  for (i in seq_len(length(thresholds) - window_size)) {
    genes <- unique(unlist(strsplit(rownames(ranked_features[i:(i + window_size), , drop = FALSE]), "→")))
    gene_entrez <- mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")
    enrichment <- enrichDGN(unname(gene_entrez), universe = unname(rbp_background))
    result <- enrichment@result
    
    if (nrow(result) == 0) {
      window_pvalues[[i]] <- 0
    } else {
      liver_hits <- unique(c(grep("liver", result$Description), grep("Liver", result$Description), grep("hepato", result$Description), grep("Hepato", result$Description)))
      window_pvalues[[i]] <- result$p.adjust[liver_hits]
      names(window_pvalues[[i]]) <- result$Description[liver_hits]
    }
  }
  
  p_top <- pval_mean(top_pvalues)
  p_window <- pval_mean(window_pvalues)
  
  p1 <- ggplot(data.frame(p = -log10(p_top), threshold = thresholds[seq_along(p_top)]), aes(x = threshold, y = p)) +
    geom_bar(stat = "identity", fill = "#7A7A7A") +
    theme_minimal() +
    theme(legend.position = "none", text = element_text(size = 15)) +
    xlab("number of links") +
    ylab(expression(-log[10]("P value"))) +
    stat_peaks(col = "#0F61AF", span = 30, geom = "text_s", ignore_threshold = 0.05, size = 6, point.padding = 0.7)
  
  p2 <- ggplot(data.frame(p = p_window, threshold = thresholds[seq_along(p_window)]), aes(x = threshold, y = p)) +
    geom_bar(stat = "identity", fill = "#7EAFC7") +
    geom_hline(yintercept = DISEASE_PVALUE_THRESHOLD, linetype = "dashed", color = "#FE3200") +
    theme_minimal() +
    xlab(paste0("sliding window of size ", window_size)) +
    ylab("P value")
  
  list(plot = ggarrange(p1, p2, nrow = 2), top_pvalues = p_top, window_pvalues = p_window)
}

# -----------------------------
# LOAD INPUTS
# -----------------------------
gtex_pheno <- fread(GTEX_PHENOTYPE_FILE)
lihc_surv <- fread(LIHC_SURVIVAL_FILE)
genes_id <- fread(GENE_PROBEMAP_FILE)
liver_expression <- readRDS(LIVER_EXPRESSION_FILE)

rbp_table <- read_excel(RBP_LIST_FILE)
cancer_lit <- read.table(CANCER_LITERATURE_FILE, sep = ",", header = TRUE)
cancer_deg <- read.table(CANCER_DEG_FILE, sep = ",", header = TRUE)
funcoup_network <- read_tsv(FUNCOUP_NETWORK_FILE, show_col_types = FALSE)
annotation <- import.gff(GENCODE_ANNOTATION_FILE)
eclip_metadata <- read_tsv(ECLIP_METADATA_FILE, show_col_types = FALSE)

# -----------------------------
# STAGE 1: EXPRESSION, SURVIVAL, AND CO-EXPRESSION VALIDATION
# -----------------------------
edges <- readRDS(ALL_INTERACTIONS_FILE)
colnames(edges) <- c("fromi", "toi", "weight", "from", "to", "color", "arrows", "arrows_type", "k562")
edges <- edges[, c("from", "to", "weight", "k562")]

analysis_results <- t(apply(edges, 1, function(row) {
  analyze_interaction_expression(row[["from"]], row[["to"]], liver_expression, lihc_surv)
}))

colnames(analysis_results) <- c("deg1", "deg2", "survg1", "survg2", "corg1", "corg2", "corpg1", "corpg2")
edges_summary <- cbind(edges, as.data.frame(analysis_results))
saveRDS(edges_summary, OUTPUT_EDGE_SUMMARY_RDS)

# -----------------------------
# STAGE 2: RBP/CANCER ANNOTATION AND FUNCOUP CLASSIFICATION
# -----------------------------
edges_annotated <- edges_summary
edges_annotated$from <- to_public_symbol(edges_annotated$from)
edges_annotated$to <- to_public_symbol(edges_annotated$to)

edges_annotated <- add_rbp_cancer_annotations(edges_annotated, rbp_table, cancer_lit, cancer_deg)

for (column in c("deg1", "deg2", "survg1", "survg2", "corpg1", "corpg2")) {
  edges_annotated[[column]] <- discretize_pvalue(edges_annotated[[column]])
}

edges_annotated$regulator_canonical <- as.numeric(as.factor(edges_annotated$regulator_canonical))
edges_annotated$target_canonical <- as.numeric(as.factor(edges_annotated$target_canonical))

edges_annotated <- add_funcoup_scores(edges_annotated, funcoup_network)
rownames(edges_annotated) <- paste0(edges_annotated$from, "-", edges_annotated$to)

funcoup_result <- impute_funcoup_classes(edges_annotated)
final_edges <- funcoup_result$edges
saveRDS(final_edges, OUTPUT_FINAL_EDGE_TABLE_RDS)

# -----------------------------
# STAGE 3: eCLIP SUPPORT
# -----------------------------
interaction_names <- do.call(rbind, strsplit(rownames(final_edges), "-")) %>% as.data.frame()
colnames(interaction_names) <- c("regulator", "target")

final_edges$eclip_reg <- mapply(
  get_eclip_pair_support,
  regulator = interaction_names$regulator,
  target = interaction_names$target,
  MoreArgs = list(metadata = eclip_metadata, annotation = annotation, eclip_dir = ECLIP_DIR, fdr_threshold = ECLIP_FDR_THRESHOLD)
)

final_edges$eclip_targ <- mapply(
  get_eclip_pair_support,
  regulator = interaction_names$target,
  target = interaction_names$regulator,
  MoreArgs = list(metadata = eclip_metadata, annotation = annotation, eclip_dir = ECLIP_DIR, fdr_threshold = ECLIP_FDR_THRESHOLD)
)

saveRDS(final_edges, OUTPUT_ECLIP_VALIDATED_RDS)

# -----------------------------
# STAGE 4: MULTI-FEATURE LINK RANKING
# -----------------------------
score_result <- build_validation_score_matrix(final_edges, N_MCLUST_GROUPS)
feature_matrix <- score_result$features
rank_score <- score_result$rank_score

ranked_features <- feature_matrix[order(rank_score, decreasing = TRUE), ]
ranked_edges <- final_edges[order(rank_score, decreasing = TRUE), ]
ranked_scores <- rank_score[order(rank_score, decreasing = TRUE)]

ranked_link_table <- as.data.frame(do.call(rbind, strsplit(rownames(ranked_features), "→")))
colnames(ranked_link_table) <- c("regulator", "target")
ranked_link_table$grn_freq <- abs(final_edges$cfreq)[order(rank_score, decreasing = TRUE)] * 8
ranked_link_table$funcoup_score <- ranked_features$fc_score
write.csv2(ranked_link_table, OUTPUT_RANKED_LINKS_CSV, row.names = FALSE)

# -----------------------------
# STAGE 5: DISEASE ENRICHMENT ACROSS RANKED LINKS
# -----------------------------
rbp_background_entrez <- mapIds(org.Hs.eg.db, rbp_table$gene_name, "ENTREZID", "SYMBOL")
rbp_background_entrez <- rbp_background_entrez[!is.na(rbp_background_entrez)]

disease_enrichment <- plot_rank_disease_enrichment(
  ranked_features = ranked_features,
  rbp_background = rbp_background_entrez,
  window_size = DISEASE_WINDOW_SIZE
)

ggsave(OUTPUT_DISEASE_ENRICHMENT_FIGURE, plot = disease_enrichment$plot, width = 10, height = 10, bg = "white")

top_indices <- seq_len(TOP_N_INTERACTIONS)
top_interactions <- as.data.frame(do.call(rbind, strsplit(rownames(ranked_features[top_indices, ]), "→")))
colnames(top_interactions) <- c("regulator", "target")
write.csv2(top_interactions, OUTPUT_TOP_INTERACTIONS_CSV, row.names = FALSE)

# -----------------------------
# STAGE 6: VALIDATION FEATURE TABLE PLOT
# -----------------------------
top_features <- ranked_features[top_indices, ]
top_edges <- ranked_edges[top_indices, ]
top_features[top_features == 0] <- NaN

top_features <- top_features[, c(order(colMeans(top_features[, -ncol(top_features)], na.rm = TRUE), decreasing = TRUE), ncol(top_features))]
top_features <- top_features[, c(order(apply(top_features[, -ncol(top_features)], 2, function(x) sum(is.na(x)))), ncol(top_features))]

label_df <- data.frame(
  GRN_frequency = format(abs(top_edges$cfreq) * 8),
  in_k562 = gsub("2", "c", gsub("1", "unc", top_edges$k562i)),
  DEG_LIHC = paste0(pvalue_score_to_symbol(top_edges$regde), "→", pvalue_score_to_symbol(top_edges$tarde)),
  alter_survival = paste0(pvalue_score_to_symbol(top_edges$regsur), "→", pvalue_score_to_symbol(top_edges$tarsurv)),
  coexpression_change = paste0("C:", round(top_edges$corLIHC, 2), ",H:", round(top_edges$corGTEx, 2)),
  times_listed_as_RBP = gsub("NA", "Ø", paste0(round(top_edges$regtlrbp, 2), "→", round(top_edges$tartlrbp, 2))),
  in_literature_as_cancer = gsub("NA", "Ø", paste0(round(top_edges$regclit, 2), "→", round(top_edges$tarclit, 2))),
  in_cancers_as_DEG = gsub("NA", "Ø", paste0(round(top_edges$regcdeg, 2), "→", round(top_edges$tarcdeg, 2))),
  RBD = gsub("NA", "Ø", paste0(gsub("2", "nc", gsub("1", "c", top_edges$regcan)), "→", gsub("2", "nc", gsub("1", "c", top_edges$tarcan)))),
  regulator_eCLIP = format(round(top_edges$eclip_reg, 1)),
  target_eCLIP = format(round(top_edges$eclip_targ, 1)),
  fc_score = top_edges$fc_score_binned
)
rownames(label_df) <- rownames(top_features)

plot_data <- top_features[, -ncol(top_features)] %>%
  rownames_to_column("rowname") %>%
  gather(colname, value, -rowname)

label_data <- label_df[, colnames(top_features)[-ncol(top_features)], drop = FALSE] %>%
  rownames_to_column("rowname") %>%
  gather(colname, label, -rowname)

plot_data$label <- label_data$label
plot_data$colname <- factor(plot_data$colname, levels = colnames(top_features))
plot_data$rowname <- factor(plot_data$rowname, levels = rev(unique(plot_data$rowname)))

text_colors <- ifelse(is.na(plot_data$value), "grey30", "white")
text_colors[plot_data$label == "NaN"] <- "white"

p_heatmap <- ggplot(plot_data, aes(x = colname, y = rowname, fill = value)) +
  geom_tile(colour = "#636363", width = 1) +
  geom_text(aes(label = label), color = text_colors, size = 3) +
  annotate(geom = "text", x = length(unique(plot_data$colname)) + 0.6, label = rev(top_features$fc_score), y = seq_len(nrow(top_features)), hjust = 0, col = "grey30") +
  scale_fill_gradient(high = "#ff6200", low = "#2275EC", na.value = "white") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    legend.position = "left",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  coord_cartesian(xlim = c(0, length(unique(plot_data$colname)) + 4), clip = "off") +
  xlab("validation features") +
  ylab("")

rank_bar_data <- data.frame(
  val = ranked_scores[top_indices],
  inters = factor(rownames(ranked_features[top_indices, ]), levels = rownames(ranked_features[top_indices, ]))
)

p_rankbar <- ggplot(rank_bar_data, aes(x = rev(val), y = inters)) +
  geom_bar(stat = "identity", colour = "white", fill = "#636363") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylab("") +
  xlab("average score")

validation_plot <- ggarrange(p_heatmap, p_rankbar, ncol = 2, labels = c("A", ""), widths = c(1.2, 0.2))
ggsave(OUTPUT_VALIDATION_HEATMAP, plot = validation_plot, width = 15, height = 5, dpi = 700, bg = "white")

# -----------------------------
# STAGE 7: NETWORK VISUALIZATION BY FUNCOUP CLASS
# -----------------------------
network_df <- as.data.frame(do.call(rbind, strsplit(rownames(top_features), "→")))
colnames(network_df) <- c("from", "to")
network_df$fcscorelabel <- top_features$fc_score

p_low <- ggplot(
  ggnetwork(network_df[network_df$fcscorelabel %in% c("low", "low*"), ], arrow.gap = 0.015),
  aes(x, y, xend = xend, yend = yend)
) +
  geom_edges(aes(color = fcscorelabel), arrow = arrow(length = unit(4, "pt"), type = "closed"), curvature = 0.3, ncp = 10, linewidth = 1) +
  scale_color_manual(values = c("#07756b", "#91e3db"), name = "FunCoup score") +
  geom_nodes(size = 6, color = "#999999") +
  geom_nodetext(aes(label = vertex.names), alpha = 0.5, fontface = "bold", color = "black", size = 4, nudge_y = 0.03) +
  theme_blank(base_size = 12)

p_high <- ggplot(
  ggnetwork(network_df[network_df$fcscorelabel %in% c("high", "high*"), ], arrow.gap = 0.015),
  aes(x, y, xend = xend, yend = yend)
) +
  geom_edges(aes(color = fcscorelabel), arrow = arrow(length = unit(4, "pt"), type = "closed"), curvature = 0.3, ncp = 10, linewidth = 1) +
  scale_color_manual(values = c("#8c1aa1", "#f5b8ff"), name = "FunCoup score") +
  geom_nodes(size = 6, color = "#999999") +
  geom_nodetext(aes(label = vertex.names), alpha = 0.5, fontface = "bold", color = "black", size = 4, nudge_y = 0.03) +
  theme_blank(base_size = 12)

network_plot <- ggarrange(p_high, p_low, nrow = 2, labels = c("A", "B"))
ggsave(OUTPUT_FUNCOUP_NETWORK_FIGURE, plot = network_plot, width = 15, height = 15, dpi = 700, bg = "white")

# ============================================================
# END OF SCRIPT
# ============================================================
