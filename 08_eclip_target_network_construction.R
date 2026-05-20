# ============================================================
# Script: 08_eclip_target_network_construction.R
# Purpose: Construct an eCLIP-derived regulator-target adjacency
#          matrix by mapping significant eCLIP peaks to protein-coding
#          genes and retaining targets that are also assayed regulators.
# ============================================================

# -----------------------------
# CONFIGURATION
# -----------------------------
INPUT_DIR <- "path/to/input/"
OUTPUT_DIR <- "path/to/output/"

# Input files and directories
ECLIP_DIR <- file.path(INPUT_DIR, "eCLIP_data/")
ECLIP_METADATA_FILE <- file.path(ECLIP_DIR, "metadata.tsv")
GENCODE_ANNOTATION_FILE <- file.path(INPUT_DIR, "gencode.v45.annotation.gtf.gz")

# Output files
OUTPUT_TARGET_LIST_RDS <- file.path(OUTPUT_DIR, "eclip_target_lists.rds")
OUTPUT_EDGE_LIST_CSV <- file.path(OUTPUT_DIR, "eclip_regulator_target_edges.csv")
OUTPUT_ADJACENCY_MATRIX_CSV <- file.path(OUTPUT_DIR, "eclip_adjacency_matrix.csv")

# Analysis parameters
FDR_THRESHOLD <- 0.01
TARGET_FILE_FORMAT <- "bigBed narrowPeak"

# -----------------------------
# LIBRARIES
# -----------------------------
library(tidyverse)
library(igraph)
library(rtracklayer)
library(data.table)
library(RCAS)

# -----------------------------
# FUNCTIONS
# -----------------------------

# Load eCLIP metadata and genome annotation.
load_eclip_inputs <- function() {
  list(
    metadata = read_tsv(ECLIP_METADATA_FILE, show_col_types = FALSE),
    annotation = import.gff(GENCODE_ANNOTATION_FILE)
  )
}

# Extract the clean target name used in ENCODE metadata.
clean_experiment_target <- function(x) {
  gsub("-human", "", x)
}

# Extract significant protein-coding gene targets for one eCLIP experiment target.
get_eclip_targets <- function(regulator, metadata, annotation, eclip_dir, fdr_threshold = 0.01) {
  metadata_subset <- metadata[clean_experiment_target(metadata$`Experiment target`) == regulator, ]
  
  if (nrow(metadata_subset) == 0) {
    return(character())
  }
  
  file_ids <- metadata_subset$`File accession`[
    metadata_subset$`File format` == TARGET_FILE_FORMAT
  ]
  
  if (length(file_ids) == 0) {
    return(character())
  }
  
  all_targets <- list()
  
  for (file_id in file_ids) {
    peak_file <- file.path(eclip_dir, paste0(file_id, ".bigBed"))
    
    if (!file.exists(peak_file)) {
      warning("Missing eCLIP bigBed file: ", peak_file)
      next
    }
    
    peak_data <- import.bb(peak_file)
    names(mcols(peak_data)) <- c("name", "score", "signalValue", "pValue", "qValue", "peak")
    
    overlaps <- as.data.table(queryGff(queryRegions = peak_data, gffData = annotation))
    
    overlaps <- overlaps[
      gene_type == "protein_coding" &
        type == "gene"
    ]
    
    overlaps <- overlaps[
      p.adjust(10^(-as.numeric(query_pValue)), method = "fdr") <= fdr_threshold
    ]
    
    all_targets[[file_id]] <- overlaps$gene_name
  }
  
  unique(unlist(all_targets))
}

# Build target lists for all regulators represented in the metadata.
build_eclip_target_lists <- function(metadata, annotation, eclip_dir, fdr_threshold = 0.01) {
  regulators <- clean_experiment_target(metadata$`Experiment target`) %>% unique()
  
  target_lists <- lapply(regulators, function(regulator) {
    targets <- get_eclip_targets(
      regulator = regulator,
      metadata = metadata,
      annotation = annotation,
      eclip_dir = eclip_dir,
      fdr_threshold = fdr_threshold
    )
    
    # Keep only targets that are also present as assayed regulators.
    intersect(targets, regulators)
  })
  
  names(target_lists) <- regulators
  target_lists
}

# Convert target lists to an edge list.
target_lists_to_edges <- function(target_lists) {
  edge_df <- imap_dfr(target_lists, function(targets, regulator) {
    if (length(targets) == 0) {
      return(data.frame(regulator = character(), target = character()))
    }
    
    data.frame(
      regulator = rep(regulator, length(targets)),
      target = targets,
      stringsAsFactors = FALSE
    )
  })
  
  edge_df %>% distinct()
}

# Build adjacency matrix from edge list.
build_adjacency_matrix <- function(edge_df) {
  graph_obj <- graph_from_data_frame(edge_df, directed = TRUE)
  as_adjacency_matrix(graph_obj, sparse = FALSE) %>% as.matrix()
}

# -----------------------------
# RUN ANALYSIS
# -----------------------------
inputs <- load_eclip_inputs()

eclip_target_lists <- build_eclip_target_lists(
  metadata = inputs$metadata,
  annotation = inputs$annotation,
  eclip_dir = ECLIP_DIR,
  fdr_threshold = FDR_THRESHOLD
)

edge_list <- target_lists_to_edges(eclip_target_lists)
adjacency_matrix <- build_adjacency_matrix(edge_list)

# Summary statistic: mean number of retained targets per regulator.
mean_targets_per_regulator <- mean(lengths(eclip_target_lists))
message("Mean retained eCLIP targets per regulator: ", round(mean_targets_per_regulator, 3))

# -----------------------------
# SAVE OUTPUTS
# -----------------------------
saveRDS(eclip_target_lists, OUTPUT_TARGET_LIST_RDS)
write.csv(edge_list, OUTPUT_EDGE_LIST_CSV, row.names = FALSE)
write.csv(adjacency_matrix, OUTPUT_ADJACENCY_MATRIX_CSV, row.names = TRUE)

# ============================================================
# END OF SCRIPT
# ============================================================