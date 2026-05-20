# ------------------------------------------------------------
# Calculate overlap statistics between two gene sets.
#
# Computes:
#   - observed overlap (intersection size),
#   - expected overlap under random sampling,
#   - hypergeometric p-value for overlap significance.
#
# Parameters:
#   list1         : First input vector (e.g., gene set A)
#   list2         : Second input vector (e.g., gene set B)
#   total_size    : Size of the background universe
#   lower_tail    : If TRUE, computes P[X <= x];
#                   if FALSE, computes P[X > x]
#   adjust_upper  : If TRUE and lower_tail = FALSE,
#                   adjusts to compute P[X >= x]
#
# Returns:
#   A named list containing:
#     - actual_overlap
#     - expected_overlap
#     - p_value
# ------------------------------------------------------------

calculate_overlap_and_pvalue <- function(
    list1,
    list2,
    total_size,
    lower_tail = TRUE,
    adjust_upper = FALSE
) {
  
  # Validate input.
  if (total_size <= 0) {
    stop("total_size must be greater than zero.")
  }
  
  if (length(list1) > total_size || length(list2) > total_size) {
    stop("Input list sizes cannot exceed total_size.")
  }
  
  # Calculate observed overlap.
  actual_overlap <- length(intersect(list1, list2))
  
  # Calculate expected overlap under independence.
  # Cast to numeric to avoid integer overflow.
  expected_overlap <- as.numeric(length(list1)) *
    length(list2) / total_size
  
  # Optional adjustment for upper-tail probability:
  # convert P[X > x] to P[X >= x].
  adjustment_value <- 0
  
  if (adjust_upper && !lower_tail) {
    adjustment_value <- 1
    warning("Calculating upper-tail probability: P[X >= x]")
  }
  
  # Hypergeometric test.
  overlap_p_value <- phyper(
    q = actual_overlap - adjustment_value,
    m = length(list1),
    n = total_size - length(list1),
    k = length(list2),
    lower.tail = lower_tail
  )
  
  # Return results.
  list(
    actual_overlap = actual_overlap,
    expected_overlap = expected_overlap,
    p_value = overlap_p_value
  )
}