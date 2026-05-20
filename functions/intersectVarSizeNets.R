# ------------------------------------------------------------
# Compare links between two networks with different gene sets.
#
# For each non-zero link in network_1, this function checks whether
# the same regulator-target pair is available and present in network_2.
#
# Assumptions:
#   - Both networks are adjacency matrices
#   - Row and column names contain gene names
#   - Row/column names are internally consistent within each network
#
# Output values:
#   common_link      : link exists in both networks
#   uncommon_link    : link exists in network_1 but not in network_2
#   unavailable_link : comparison is impossible because one or both genes
#                      are absent from network_2
#
# Returns:
#   A matrix with the same dimensions as network_1.
# ------------------------------------------------------------

intersect_variable_size_networks <- function(network_1, network_2) {
  
  # Validate input.
  if (is.null(rownames(network_1)) || is.null(colnames(network_1))) {
    stop("network_1 must have row and column names.")
  }
  
  if (is.null(rownames(network_2)) || is.null(colnames(network_2))) {
    stop("network_2 must have row and column names.")
  }
  
  if (!identical(rownames(network_1), colnames(network_1))) {
    warning("network_1 row and column names are not identical.")
  }
  
  if (!identical(rownames(network_2), colnames(network_2))) {
    warning("network_2 row and column names are not identical.")
  }
  
  genes_1 <- colnames(network_1)
  genes_2 <- colnames(network_2)
  
  # Output matrix keeps original network_1 structure.
  comparison_matrix <- matrix(
    NA_character_,
    nrow = nrow(network_1),
    ncol = ncol(network_1),
    dimnames = dimnames(network_1)
  )
  
  # Compare only links present in network_1.
  for (source_gene in genes_1) {
    
    targets_in_network_1 <- rownames(network_1)[
      as.logical(network_1[, source_gene])
    ]
    
    if (length(targets_in_network_1) == 0) {
      next
    }
    
    # If source gene is missing in network_2, all its links are unavailable.
    if (!source_gene %in% genes_2) {
      comparison_matrix[targets_in_network_1, source_gene] <- "unavailable_link"
      next
    }
    
    for (target_gene in targets_in_network_1) {
      
      # If target gene is missing in network_2, link cannot be compared.
      if (!target_gene %in% rownames(network_2)) {
        comparison_matrix[target_gene, source_gene] <- "unavailable_link"
        next
      }
      
      # Check whether the same link exists in network_2.
      link_value_network_2 <- network_2[target_gene, source_gene]
      
      if (is.na(link_value_network_2)) {
        comparison_matrix[target_gene, source_gene] <- "unavailable_link"
      } else if (abs(link_value_network_2) > 0) {
        comparison_matrix[target_gene, source_gene] <- "common_link"
      } else {
        comparison_matrix[target_gene, source_gene] <- "uncommon_link"
      }
    }
  }
  
  return(comparison_matrix)
}