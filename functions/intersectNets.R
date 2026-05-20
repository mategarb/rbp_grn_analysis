# ------------------------------------------------------------
# Intersect multiple inferred networks by edge frequency.
#
# Combines a list of adjacency matrices and generates a series
# of pruned consensus networks based on edge occurrence frequency.
#
# Assumptions:
#   - All input matrices have the same dimensions
#   - All matrices use the same gene ordering
#
# Parameters:
#   network_list        : List of adjacency matrices
#   add_gene_names      : Whether to assign row/column names
#   gene_names          : Vector of gene names (optional)
#   binary_input        : If TRUE, negative values are converted to 1
#                         (for undirected/binary networks)
#   binary_output       : If TRUE, output is discretized to {-1,1}
#                         instead of frequency values
#   discretize_input    : If TRUE, convert input weights to {-1,0,1}
#   remove_diagonal     : If TRUE, set diagonal values to zero
#
# Returns:
#   A named list of consensus networks at increasing
#   frequency thresholds.
#
# Example:
#   intersect_networks(
#     network_list = nets,
#     gene_names = genes
#   )
# ------------------------------------------------------------

intersect_networks <- function(
    network_list,
    add_gene_names = TRUE,
    gene_names = NULL,
    binary_input = FALSE,
    binary_output = FALSE,
    discretize_input = TRUE,
    remove_diagonal = TRUE
) {
  
  # Validate input.
  if (!is.list(network_list) || length(network_list) == 0) {
    stop("network_list must be a non-empty list of matrices.")
  }
  
  network_dimensions <- lapply(network_list, dim)
  
  if (!all(sapply(network_dimensions, identical, network_dimensions[[1]]))) {
    stop("All network matrices must have identical dimensions.")
  }
  
  # Optional discretization of input networks.
  if (discretize_input) {
    for (i in seq_along(network_list)) {
      
      if (binary_input) {
        # Convert all non-zero values to 1.
        network_list[[i]][network_list[[i]] < 0] <- 1
      } else {
        # Preserve sign information.
        network_list[[i]][network_list[[i]] < 0] <- -1
      }
      
      network_list[[i]][network_list[[i]] > 0] <- 1
    }
  }
  
  # Compute average edge frequency across all networks.
  average_network <- Reduce("+", network_list) / length(network_list)
  
  consensus_networks <- vector("list", length(network_list))
  
  # Generate thresholded consensus networks.
  for (threshold_index in seq_along(network_list)) {
    
    threshold_network <- average_network
    
    # Keep only edges present in at least threshold proportion.
    threshold_value <- threshold_index / length(network_list)
    
    threshold_network[
      abs(threshold_network) < threshold_value &
        threshold_network != 0
    ] <- 0
    
    # Optional binary output.
    if (binary_output) {
      threshold_network[threshold_network < 0] <- -1
      threshold_network[threshold_network > 0] <- 1
    }
    
    # Optional gene naming.
    if (add_gene_names) {
      if (is.null(gene_names)) {
        stop("gene_names must be provided when add_gene_names = TRUE.")
      }
      
      colnames(threshold_network) <- gene_names
      rownames(threshold_network) <- gene_names
    }
    
    # Optional diagonal cleanup.
    if (remove_diagonal) {
      diag(threshold_network) <- 0
    }
    
    consensus_networks[[threshold_index]] <- threshold_network
  }
  
  # Assign names indicating frequency thresholds.
  names(consensus_networks) <- paste0(
    seq_along(network_list),
    "/",
    length(network_list)
  )
  
  return(consensus_networks)
}