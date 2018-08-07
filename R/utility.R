
normalize_cols <- function(M, ranked = TRUE) {
  if (ranked) {
    M <- apply(M, 2, rank) - (nrow(M)+1)/2
  } else {
    M <- scale(M, scale=FALSE)
  }
  M <- scale(M, center=FALSE, apply(M, 2, function(c) sqrt(sum(c**2))))
  return(M)
}

find_subsets <- function(full_list, list_names) {
  return(sapply(list_names, function(name) full_list == name))
}

compute_aurocs <- function(votes) {
  candidate_labels <- rownames(votes)
  positives <- design_matrix(candidate_labels)
  n_positives <- colSums(positives)
  n_negatives <- nrow(positives) - n_positives
  sum_of_positive_ranks <- t(positives) %*% apply(abs(votes), MARGIN = 2, FUN = rank)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

design_matrix <- function(cell_type) {
  factors <- levels(as.factor(cell_type))
  if (length(factors) > 1) {
    result <- model.matrix(~cell_type-1)
  } else {
    result <- matrix(1, nrow = length(cell_type), ncol = 1)
  }
  colnames(result) <- factors
  return(result)
}

create_result_matrix <- function(cell_type) {
  unique_cell_type <- unique(cell_type)
  result <- matrix(0, nrow = length(unique_cell_type), ncol = length(unique_cell_type))
  rownames(result) <- unique_cell_type
  colnames(result) <- unique_cell_type
  return(result)
}

get_study_id <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "\\|"), head, 1))
}

get_cell_type <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "\\|"), tail, 1))
}
