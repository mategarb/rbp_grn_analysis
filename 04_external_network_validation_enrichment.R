# ============================================================
# Script: 04_external_network_validation_enrichment.R
# Purpose: Validate selected regulator-target interactions using
#          external liver regulatory-network resources, visualize
#          validation weights, and perform functional enrichment.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files
VALIDATION_TABLE_FILE <- file.path(INPUT_DIR, "Table_S3_validationtable_3plus.csv")
GRAND_CELL_LINE_DIR <- file.path(INPUT_DIR, "Grand_cell_lines/")
GRAND_METADATA_FILE <- file.path(INPUT_DIR, "Grand_others/meta_grand.xlsx")
TCGA_OTTER_NETWORK_FILE <- file.path(INPUT_DIR, "Grand_cell_lines/cancer_liver_otter_network.csv")
HEPG2_EXPRESSION_FILE <- file.path(INPUT_DIR, "datasets/ymatrix_hepg2.csv")
RBP_LIST_FILE <- file.path(INPUT_DIR, "rbps_info/210329_Table_S1_hRBP_list.xlsx")
MSIGDB_CANONICAL_PATHWAYS_FILE <- file.path(INPUT_DIR, "MSigDB/c2.cp.v2023.2.Hs.symbols.gmt")

# GRNdb regulon files
GRNDB_LIVER_TUMOR_FILE <- file.path(INPUT_DIR, "GRNdb/Liver-Tumor-regulons.txt")
GRNDB_LIVER_NORMAL_FILE <- file.path(INPUT_DIR, "GRNdb/Liver-Normal-regulons.txt")
GRNDB_ADULT_LIVER_FILE <- file.path(INPUT_DIR, "GRNdb/Adult-Liver-regulons.txt")
GRNDB_LIHC_TCGA_FILE <- file.path(INPUT_DIR, "GRNdb/LIHC_TCGA-regulons.txt")
GRNDB_LIVER_BLOOD_FILE <- file.path(INPUT_DIR, "GRNdb/Liver-Peripheral-blood-regulons.txt")

# Optional files used for exploratory MYC checks
MYC_TARGETS_V1_FILE <- file.path(INPUT_DIR, "MSigDB/HALLMARK_MYC_TARGETS_V1.v2023.2.Hs.tsv")
MYC_TARGETS_V2_FILE <- file.path(INPUT_DIR, "MSigDB/HALLMARK_MYC_TARGETS_V2.v2023.2.Hs.tsv")
NCRNA_GENE_SILENCING_FILE <- file.path(INPUT_DIR, "MSigDB/GOBP_REGULATORY_NCRNA_MEDIATED_GENE_SILENCING.v2023.2.Hs.tsv")

# Output files
OUTPUT_VALIDATION_FIGURE <- file.path(OUTPUT_DIR, "external_network_validation.svg")
OUTPUT_ENRICHMENT_FIGURE <- file.path(OUTPUT_DIR, "external_validation_enrichment.svg")
OUTPUT_GRAND_VALUES_RDS <- file.path(OUTPUT_DIR, "grand_validation_values.rds")
OUTPUT_GRAND_DIRECTION_RDS <- file.path(OUTPUT_DIR, "grand_validation_direction.rds")

# Analysis parameters
GRAND_N_CELL_LINES <- 24
P_VALUE_THRESHOLD <- 0.05
TOP_ENRICHMENT_TERMS <- 20

# -----------------------------
# LIBRARIES
# -----------------------------
library(data.table)
library(tidyverse)
library(readxl)
library(igraph)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
library(gprofiler2)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(stringr)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Read selected interaction pairs from the validation table.
read_interaction_pairs <- function(validation_table_file) {
  inters <- read.csv2(validation_table_file)
  inters <- do.call(rbind, strsplit(inters[, 1], "-")) %>% as.data.frame()
  colnames(inters) <- c("regulator", "target")
  inters
}

# Extract interaction weights from a long-format regulatory network.
# If the original direction is missing, the swapped direction is checked.
extract_link_weights <- function(interactions, network_long, source_col = "source", target_col = "target", value_col = "value") {
  values <- rep(NA_real_, nrow(interactions))
  direction <- rep(NA_character_, nrow(interactions))
  
  for (i in seq_len(nrow(interactions))) {
    regulator <- interactions$regulator[i]
    target <- interactions$target[i]
    
    direct_idx <- which(network_long[[source_col]] == regulator & network_long[[target_col]] == target)
    swapped_idx <- which(network_long[[source_col]] == target & network_long[[target_col]] == regulator)
    
    if (length(direct_idx) > 0) {
      values[i] <- network_long[[value_col]][direct_idx[1]]
      direction[i] <- "unchanged"
    } else if (length(swapped_idx) > 0) {
      values[i] <- network_long[[value_col]][swapped_idx[1]]
      direction[i] <- "swapped"
    }
  }
  
  list(values = values, direction = direction)
}

# Convert a matrix-style network file to long format.
read_matrix_network_long <- function(file_path, row_id_col = "Row") {
  network <- fread(file_path) %>% as.data.frame()
  melt(network, id.vars = row_id_col) %>%
    rename(source = all_of(row_id_col), target = variable, value = value)
}

# Validate all selected links across GRAND cell-line networks.
validate_against_grand_cell_lines <- function(interactions, grand_dir, metadata_file, n_cell_lines = 24) {
  files <- list.files(path = grand_dir)
  selected_files <- files[seq_len(min(n_cell_lines, length(files)))]
  
  value_list <- list()
  direction_list <- list()
  
  for (j in seq_along(selected_files)) {
    network_long <- read_matrix_network_long(file.path(grand_dir, selected_files[j]), row_id_col = "Row")
    extracted <- extract_link_weights(interactions, network_long)
    value_list[[j]] <- extracted$values
    direction_list[[j]] <- extracted$direction
  }
  
  metadata <- read_xlsx(metadata_file, col_names = FALSE) %>% as.matrix()
  
  value_df <- do.call(rbind.data.frame, value_list) %>% t()
  rownames(value_df) <- paste0(interactions$regulator, "-", interactions$target)
  colnames(value_df) <- metadata[match(gsub(".csv", "", selected_files), metadata[, 1]), 2]
  
  value_df <- value_df[-which(apply(value_df, 1, function(x) all(is.na(x)))), , drop = FALSE]
  
  annotation_df <- data.frame(
    Sex = metadata[match(gsub(".csv", "", selected_files), metadata[, 1]), 4],
    Ethnicity = metadata[match(gsub(".csv", "", selected_files), metadata[, 1]), 5],
    Type = metadata[match(gsub(".csv", "", selected_files), metadata[, 1]), 8]
  )
  rownames(annotation_df) <- colnames(value_df)
  
  list(values = value_df, directions = direction_list, annotation = annotation_df)
}

# Plot GRAND validation heatmap.
plot_grand_heatmap <- function(value_df, annotation_df) {
  annotation_colors <- list(
    Type = c("Metastasis" = "#7d5f16", "Primary" = "#c99e38", "-" = "#F1F1F1"),
    Sex = c("Female" = "#06bf3a", "Male" = "#f2f542"),
    Ethnicity = c("african_american" = "#a1a3e3", "asian" = "#62638a", "caucasian" = "#242433", "-" = "#F1F1F1")
  )
  
  as.ggplot(
    pheatmap(
      value_df,
      scale = "none",
      annotation_col = annotation_df,
      annotation_colors = annotation_colors,
      color = colorRampPalette(c("#ad42f5", "white", "#e35e1b"))(12)
    )
  )
}

# Create a bar plot of validated link weights.
plot_link_weights <- function(weight_df, title = NULL, gradient_low = "#e8baa2", gradient_high = "#e35e1b") {
  ggplot(weight_df, aes(x = reorder(link, -weight), y = weight, fill = weight)) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_gradient(low = gradient_low, high = gradient_high) +
    theme_classic() +
    theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("link") +
    ggtitle(title)
}

# Validate links against a GRNdb-style edge table.
validate_against_edge_table <- function(interactions, edge_table, source_col, target_col, value_col) {
  network_long <- edge_table[, c(source_col, target_col, value_col)]
  colnames(network_long) <- c("source", "target", "value")
  
  extracted <- extract_link_weights(interactions, network_long)
  names(extracted$values) <- paste0(interactions$regulator, "-", interactions$target)
  
  values <- extracted$values[!is.na(extracted$values)]
  data.frame(link = names(values), weight = unname(values))
}

# Format enrichment results for plotting.
format_enrichment_result <- function(enrich_result, p_threshold = 0.05, top_n = 20) {
  result <- enrich_result@result
  result <- result[result$p.adjust <= p_threshold, ]
  
  if (nrow(result) == 0) return(data.frame())
  
  data.frame(
    pvalue = result$p.adjust,
    term = result$Description,
    count = as.numeric(sapply(strsplit(result$GeneRatio, "/"), `[[`, 1))
  ) %>%
    head(top_n)
}

# Plot enrichment bar chart.
plot_enrichment_terms <- function(enrichment_df, title) {
  if (nrow(enrichment_df) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste(title, "- no significant terms")))
  }
  
  enrichment_df$term <- factor(enrichment_df$term, levels = unique(enrichment_df$term))
  
  ggplot(enrichment_df, aes(x = term, y = count)) +
    geom_bar(aes(fill = pvalue), stat = "identity") +
    scale_fill_gradient(low = "#f35c87", high = "#5CB8F3", na.value = NA, name = "P value") +
    coord_flip() +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 50), limits = rev) +
    xlab("") +
    theme_classic() +
    theme(text = element_text(size = 16)) +
    ggtitle(title)
}

# -----------------------------
# LOAD SELECTED INTERACTIONS
# -----------------------------
interactions <- read_interaction_pairs(VALIDATION_TABLE_FILE)

# -----------------------------
# GRAND CELL-LINE VALIDATION
# -----------------------------
grand_validation <- validate_against_grand_cell_lines(
  interactions = interactions,
  grand_dir = GRAND_CELL_LINE_DIR,
  metadata_file = GRAND_METADATA_FILE,
  n_cell_lines = GRAND_N_CELL_LINES
)

saveRDS(grand_validation$values, OUTPUT_GRAND_VALUES_RDS)
saveRDS(grand_validation$directions, OUTPUT_GRAND_DIRECTION_RDS)

p_grand <- plot_grand_heatmap(
  value_df = grand_validation$values,
  annotation_df = grand_validation$annotation
)

# -----------------------------
# TCGA/OTTER VALIDATION
# -----------------------------
tcga_network <- fread(TCGA_OTTER_NETWORK_FILE) %>% as.data.frame()

# Convert gene symbols to Entrez IDs when needed, matching the original workflow.
gene_symbols <- colnames(tcga_network)
gene_entrez <- gconvert(gene_symbols, organism = "hsapiens", target = "ENTREZGENE", filter_na = FALSE)$target
colnames(tcga_network) <- gene_entrez[match(gene_symbols, colnames(tcga_network))]

# Keep HepG2 expression loading for reproducibility and symbol compatibility checks.
hepg2_expression <- read.csv2(HEPG2_EXPRESSION_FILE, header = TRUE, sep = "\t", row.names = 1)
hepg2_expression <- as.data.frame(sapply(hepg2_expression, as.numeric))

network_long <- melt(tcga_network, id.vars = "V1") %>%
  rename(source = V1, target = variable, value = value)

tcga_extracted <- extract_link_weights(interactions, network_long)
names(tcga_extracted$values) <- paste0(interactions$regulator, "-", interactions$target)

tcga_values <- tcga_extracted$values[!is.na(tcga_extracted$values)]
tcga_values <- tcga_values[tcga_values != 0]
tcga_weight_df <- data.frame(link = names(tcga_values), weight = unname(tcga_values))

p_tcga <- plot_link_weights(
  tcga_weight_df,
  title = "TCGA/OTTER liver network",
  gradient_low = "#ad42f5",
  gradient_high = "#e35e1b"
)

# -----------------------------
# GRNdb VALIDATION
# -----------------------------
grndb_liver_tumor <- read.table(GRNDB_LIVER_TUMOR_FILE, sep = "\t", header = TRUE)
grndb_liver_normal <- read.table(GRNDB_LIVER_NORMAL_FILE, sep = "\t", header = TRUE)
grndb_adult_liver <- read.table(GRNDB_ADULT_LIVER_FILE, sep = "\t", header = TRUE)
grndb_lihc_tcga <- read.table(GRNDB_LIHC_TCGA_FILE, sep = "\t", header = TRUE)
grndb_liver_blood <- read.table(GRNDB_LIVER_BLOOD_FILE, sep = "\t", header = TRUE)

# GRNdb columns used in the original script: regulator, target, weight.
tumor_weight_df <- validate_against_edge_table(interactions, grndb_liver_tumor, 1, 2, 5)
normal_weight_df <- validate_against_edge_table(interactions, rbind(grndb_liver_normal[, c(1, 2, 5)], grndb_adult_liver[, c(1, 2, 5)]), 1, 2, 3)
blood_weight_df <- validate_against_edge_table(interactions, grndb_liver_blood, 1, 2, 5)

p_tumor <- plot_link_weights(tumor_weight_df, title = "GRNdb liver tumor")
p_normal <- plot_link_weights(normal_weight_df, title = "GRNdb liver normal")
p_blood <- plot_link_weights(blood_weight_df, title = "GRNdb liver peripheral blood")

# Combine validation plots.
p_top <- ggarrange(p_grand, p_tcga, ncol = 2, labels = "AUTO", widths = c(1, 0.75))
p_bottom <- ggarrange(p_tumor, p_normal, p_blood, ncol = 3, labels = c("C", "D", "E"))
p_validation <- ggarrange(p_top, p_bottom, nrow = 2)

ggsave(OUTPUT_VALIDATION_FIGURE, plot = p_validation, width = 14, height = 12)

# -----------------------------
# FUNCTIONAL ENRICHMENT OF VALIDATED GENES
# -----------------------------
rbp_table <- read_excel(RBP_LIST_FILE)
validated_genes <- strsplit(rownames(grand_validation$values), "-") %>% unlist() %>% unique()
validated_rbps <- intersect(rbp_table$gene_name, validated_genes)

# GO Biological Process enrichment.
validated_ensembl <- mapIds(org.Hs.eg.db, validated_rbps, "ENSEMBL", "SYMBOL")

go_bp <- enrichGO(
  gene = unname(validated_ensembl),
  universe = rbp_table$GRCh38.p7_ensembl_ID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  readable = TRUE,
  pAdjustMethod = "fdr",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

p_go <- plot_enrichment_terms(
  format_enrichment_result(go_bp, P_VALUE_THRESHOLD, TOP_ENRICHMENT_TERMS),
  title = "GO Biological Process"
)

# DisGeNET enrichment.
validated_entrez <- mapIds(org.Hs.eg.db, validated_rbps, "ENTREZID", "SYMBOL")
background_entrez <- mapIds(org.Hs.eg.db, rbp_table$gene_name, "ENTREZID", "SYMBOL")
validated_entrez <- validated_entrez[!is.na(validated_entrez)]
background_entrez <- background_entrez[!is.na(background_entrez)]

dgn <- enrichDGN(
  gene = unname(validated_entrez),
  universe = unname(background_entrez),
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

p_dgn <- plot_enrichment_terms(
  format_enrichment_result(dgn, P_VALUE_THRESHOLD, TOP_ENRICHMENT_TERMS),
  title = "DisGeNET"
)

# Canonical pathway enrichment.
canonical_pathways <- read.gmt(MSIGDB_CANONICAL_PATHWAYS_FILE)
canonical_enrich <- enricher(
  gene = validated_rbps,
  universe = rbp_table$gene_name,
  TERM2GENE = canonical_pathways,
  pAdjustMethod = "fdr"
)

pathway_df <- format_enrichment_result(canonical_enrich, P_VALUE_THRESHOLD, TOP_ENRICHMENT_TERMS)
if (nrow(pathway_df) > 0) {
  pathway_terms <- gsub("_", " ", pathway_df$term) %>% str_to_title()
  pathway_df$term <- str_remove(pathway_terms, "(\\w+\\s+){1}")
}

p_pathways <- plot_enrichment_terms(pathway_df, title = "Canonical Pathways")

p_enrichment <- ggarrange(
  p_pathways + theme(text = element_text(size = 18)),
  p_dgn + theme(text = element_text(size = 18)),
  nrow = 2,
  labels = "AUTO"
)

ggsave(OUTPUT_ENRICHMENT_FIGURE, plot = p_enrichment, width = 10, height = 14)

# -----------------------------
# OPTIONAL CYTOSCAPE EXPORT DATA
# -----------------------------
# These tables can be imported into Cytoscape manually or used with RCy3.
cytoscape_edges <- data.frame(
  source = sapply(strsplit(rownames(grand_validation$values), "-"), `[[`, 1),
  target = sapply(strsplit(rownames(grand_validation$values), "-"), `[[`, 2),
  weight = rowMeans(grand_validation$values, na.rm = TRUE)
)

cytoscape_edges$interaction <- ifelse(cytoscape_edges$weight < 0, "inhibition", "activation")
cytoscape_edges$relative_weight <- round(abs(cytoscape_edges$weight) / mean(abs(cytoscape_edges$weight), na.rm = TRUE), digits = 1)

cytoscape_nodes <- data.frame(
  id = unique(c(cytoscape_edges$source, cytoscape_edges$target)),
  group = "node",
  stringsAsFactors = FALSE
)

# -----------------------------
# OPTIONAL MYC-RELATED GENE-SET CHECKS
# -----------------------------
myc_rbps_set_1 <- c("PKM", "SERBP1", "EIF4G1", "XPO1", "FASTKD1", "RACK1", "NPM1", "HNRNPK", "YWHAG", "HSPD1", "SUPV3L1")
myc_rbps_set_2 <- c("HNRNPLL", "SF3B4", "AKAP8L", "NPM1", "RPL23A", "ADAR", "RACK1", "SUPV3L1", "SRSF3", "XRN1", "EIF3D", "SUB1", "EIF2S2", "IGF2BP1", "PCBP2", "SRSF1", "FAM120A", "CCAR1")

if (file.exists(MYC_TARGETS_V1_FILE) && file.exists(MYC_TARGETS_V2_FILE)) {
  myc_v1 <- readr::read_tsv(MYC_TARGETS_V1_FILE, show_col_types = FALSE)
  myc_v2 <- readr::read_tsv(MYC_TARGETS_V2_FILE, show_col_types = FALSE)
  
  myc_v1_genes <- myc_v1[17, 2] %>% as.matrix() %>% as.character() %>% strsplit(",") %>% unlist()
  myc_v2_genes <- myc_v2[17, 2] %>% as.matrix() %>% as.character() %>% strsplit(",") %>% unlist()
  
  myc_overlap_v1 <- intersect(unique(c(myc_rbps_set_1, myc_rbps_set_2)), myc_v1_genes)
  myc_overlap_v2 <- intersect(unique(c(myc_rbps_set_1, myc_rbps_set_2)), myc_v2_genes)
}

if (file.exists(NCRNA_GENE_SILENCING_FILE)) {
  ncrna_silencing <- readr::read_tsv(NCRNA_GENE_SILENCING_FILE, show_col_types = FALSE)
  ncrna_silencing_genes <- ncrna_silencing[17, 2] %>% as.matrix() %>% as.character() %>% strsplit(",") %>% unlist()
  ncrna_silencing_overlap <- intersect(c("PABPC1", "UTP3", "TIAL1", "LIN28B", "DDX21", "PPIL4"), ncrna_silencing_genes)
}

# ============================================================
# END OF SCRIPT
# ============================================================