# ============================================================
# Script: 02_shared_target_overlap_analysis.R
# Purpose: Analyze shared targets between regulator pairs using
#          eCLIP and RAP-seq datasets and perform pathway enrichment.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
ECLIP_DIR <- file.path(INPUT_DIR, "eclip_data/")
RAPSEQ_DIR <- file.path(INPUT_DIR, "rapseq_peaks/")
ANNOTATION_FILE <- file.path(INPUT_DIR, "gencode_annotation.gtf.gz")
METADATA_FILE <- file.path(INPUT_DIR, "metadata.tsv")
SUPP_TABLE <- file.path(INPUT_DIR, "Table_S3_validationtable.csv")
MSIGDB_FILE <- file.path(INPUT_DIR, "c2.cp.v2023.2.Hs.symbols.gmt")
OUTPUT_DIR <- "path/to/output/"

# -----------------------------
# LIBRARIES
# -----------------------------
library(tidyverse)
library(data.table)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggvenn)
library(ggpubr)
library(readxl)
library(rtracklayer)
library(GenomicRanges)
library(RCAS)
library(clusterProfiler)
library(qusage)

# -----------------------------
# INPUT DATA
# -----------------------------
# Validation interaction table
inters <- read.csv2(SUPP_TABLE)

# eCLIP metadata
meta <- read_tsv(METADATA_FILE)

# RBP annotation table
rbp_table <- read_excel(file.path(INPUT_DIR, "RBP_list.xlsx"))

# Genome annotation
annots <- import.gff(ANNOTATION_FILE)

# -----------------------------
# PREPARE INTERACTION NETWORK
# -----------------------------
interaction_pairs <- do.call(rbind, strsplit(inters[,1], "-"))
interaction_graph <- graph.edgelist(interaction_pairs, directed = TRUE)
interaction_matrix <- as_adjacency_matrix(interaction_graph) %>% as.matrix()

# Keep only experimentally supported interactions
validated_inters <- inters[
  rowSums(inters[, c("eclip_reg", "eclip_targ", "rapseq_reg", "rapseq_targ")], na.rm = TRUE) != 0,
]

regulators <- strsplit(validated_inters[,1], "-") %>% unlist() %>% unique()

# -----------------------------
# COLLECT eCLIP TARGETS
# -----------------------------
get_eclip_targets <- function(gene_name, meta, annots, eclip_dir) {
  meta_sub <- meta[gsub("-human", "", meta$`Experiment target`) == gene_name, ]
  
  if (nrow(meta_sub) == 0) return(character())
  
  file_id <- meta_sub$`File accession`[
    meta_sub$`File format` == "bigBed narrowPeak"
  ][1]
  
  peak_data <- import.bb(file.path(eclip_dir, paste0(file_id, ".bigBed")))
  names(mcols(peak_data)) <- c("name", "score", "signalValue", "pValue", "qValue", "peak")
  
  overlaps <- as.data.table(queryGff(queryRegions = peak_data, gffData = annots))
  
  overlaps <- overlaps[
    gene_type == "protein_coding" &
      type == "gene"
  ]
  
  overlaps <- overlaps[
    p.adjust(10^(-as.numeric(query_pValue)), method = "fdr") <= 0.05
  ]
  
  unique(overlaps$gene_name)
}

eclip_targets <- lapply(regulators, get_eclip_targets,
                        meta = meta,
                        annots = annots,
                        eclip_dir = ECLIP_DIR)
names(eclip_targets) <- regulators

# -----------------------------
# COLLECT RAP-seq TARGETS
# -----------------------------
rap_files <- list.files(RAPSEQ_DIR)
rap_regulators <- gsub("\\..*", "", rap_files)
rap_regulators <- intersect(rap_regulators, regulators)

get_rapseq_targets <- function(gene_name, rapseq_dir) {
  file <- file.path(rapseq_dir, paste0(gene_name, ".peaks.txt"))
  peaks <- read.table(file, header = TRUE)
  peaks <- peaks[peaks$gene_type == "protein_coding", ]
  unique(peaks$gene_name)
}

rapseq_targets <- lapply(rap_regulators, get_rapseq_targets,
                         rapseq_dir = RAPSEQ_DIR)
names(rapseq_targets) <- rap_regulators

# -----------------------------
# OVERLAP STATISTICS
# -----------------------------
edge_pairs <- do.call(rbind, strsplit(validated_inters[,1], "-"))
background_total <- length(unique(c(unlist(eclip_targets), unlist(rapseq_targets))))

compute_overlap_pvalue <- function(reg1, reg2) {
  gs1 <- unique(c(eclip_targets[[reg1]], rapseq_targets[[reg1]]))
  gs2 <- unique(c(eclip_targets[[reg2]], rapseq_targets[[reg2]]))
  
  if (length(gs1) == 0 || length(gs2) == 0) {
    return(c(NA, NA, NA, NA))
  }
  
  overlap_n <- length(intersect(gs1, gs2))
  
  contingency <- matrix(c(
    overlap_n,
    length(gs1) - overlap_n,
    length(gs2) - overlap_n,
    background_total - length(gs1) - length(gs2) + overlap_n
  ), nrow = 2)
  
  pval <- fisher.test(contingency, alternative = "greater")$p.value
  
  c(pval, overlap_n, length(gs1), length(gs2))
}

results <- apply(edge_pairs, 1, function(x) compute_overlap_pvalue(x[1], x[2]))
results <- t(results)

overlap_df <- data.frame(
  regulator = edge_pairs[,1],
  target = edge_pairs[,2],
  pval = results[,1],
  overlap = results[,2],
  n_reg = results[,3],
  n_tar = results[,4]
)

overlap_df$pval[is.na(overlap_df$pval)] <- 1
overlap_df$fdr <- p.adjust(overlap_df$pval, method = "fdr")

significant_pairs <- overlap_df %>%
  filter(fdr <= 0.05) %>%
  arrange(fdr)

# -----------------------------
# VISUALIZATION
# -----------------------------
network_graph <- graph_from_data_frame(significant_pairs, directed = TRUE)
V(network_graph)$degree <- degree(network_graph)

plot_network <- ggraph(network_graph, layout = "linear") +
  geom_edge_arc(aes(width = -log10(fdr), color = -log10(fdr)), alpha = 0.8) +
  geom_node_point(aes(size = degree), color = "gray40") +
  geom_node_label(aes(label = name)) +
  theme_void()

ggsave(
  filename = file.path(OUTPUT_DIR, "shared_target_network.svg"),
  plot = plot_network,
  width = 15,
  height = 8
)

# -----------------------------
# PATHWAY ENRICHMENT
# -----------------------------
pathways <- read.gmt(MSIGDB_FILE)

run_enrichment <- function(reg1, reg2) {
  shared_genes <- intersect(
    unique(c(eclip_targets[[reg1]], rapseq_targets[[reg1]])),
    unique(c(eclip_targets[[reg2]], rapseq_targets[[reg2]]))
  )
  
  enricher(
    gene = shared_genes,
    TERM2GENE = pathways,
    pAdjustMethod = "fdr"
  )
}

# Example: top 3 pairs
top_pairs <- significant_pairs[1:3, ]
enrichment_results <- lapply(
  1:nrow(top_pairs),
  function(i) run_enrichment(top_pairs$regulator[i], top_pairs$target[i])
)

saveRDS(
  enrichment_results,
  file = file.path(OUTPUT_DIR, "shared_target_enrichment_results.rds")
)

# ============================================================
# END OF SCRIPT
# ============================================================
