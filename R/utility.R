#' Contains a collection of utility functions
#'

#' Scale matrix such that all colums sum to 0 and have l2-norm of 1
normalize_cols <- function(M, ranked = TRUE) {
  if (ranked) {
    M <- apply(M, 2, rank)
  }
  return(normalize_cols_cpp(as.matrix(M)))
}

#' Scale matrix such that all colums sum to 0 and have l2-norm of 1 (C++)
Rcpp::cppFunction('NumericMatrix normalize_cols_cpp(NumericMatrix M) {
  NumericMatrix result(M.nrow(), M.ncol());
  for (int j = 0; j < M.ncol(); j++) {
    double mean = 0;
    for (int i = 0; i < M.nrow(); i++) { mean += M(i,j); }
    mean /= M.nrow();
    for (int i = 0; i < M.nrow(); i++) { result(i,j) = M(i,j) - mean; }
    double norm = 0;
    for (int i = 0; i < M.nrow(); i++) { norm += result(i,j) * result(i,j); }
    norm = 1 / sqrt(norm);
    for (int i = 0; i < M.nrow(); i++) { result(i,j) *= norm; }
  }
  return result;
}')

#' Return binary matrix with position of elements of list_names within full_list
find_subsets <- function(full_list, list_names) {
  return(sapply(list_names, function(name) full_list == name))
}

#' Compute approximate neighbor voting without cell-cell network
#'
#' voter_id is a binary matrix indicating cell types of the voters.
#' If left empty, cell types are assumed to be the column names of the voters.
compute_votes_without_network <- function(candidates, voters, voter_id = NULL) {
  if (is.null(voter_id)) {
    voter_id <- design_matrix(colnames(voters))
  }
  raw_votes <- crossprod(candidates, voters %*% voter_id)
  return(sweep(raw_votes, 2, colSums(voter_id), FUN = "+")) /
         (c(crossprod(candidates, rowSums(voters))) + rep(ncol(voters), ncol(candidates)))
}

#' Build cell x cell correlation network from scaled matrices
build_network <- function(set_A, set_B, ranked = TRUE) {
  result <- crossprod(set_A, set_B)
  if (!ranked) { return(result) }
  A_labels <- rownames(result)
  B_labels <- colnames(result)
  result <- matrix(pseudo_rank(result), nrow = nrow(result))
  rownames(result) <- A_labels
  colnames(result) <- B_labels
  return(result)
}

#' Rank data approximately (sampling-based, calling close values as ties)
pseudo_rank <- function(x, breaks = 1000, depth = 1000) {
  m <- min(x)
  M <- max(x)
  bins <- floor((x-m) / ((1+1e-10)*(M-m)) * breaks) + 1
  if (is.null(depth)) {
    num_per_bin <- count_occurrences(bins, breaks)
    rank_per_bin <- count_to_rank(num_per_bin, length(x))
  } else {
    num_per_bin <- count_occurrences(bins[sample.int(length(x), breaks*depth)], breaks)
    rank_per_bin <- count_to_rank(num_per_bin, breaks*depth)
  }
  return(rank_per_bin[bins])
}

#' Count occurrences of integers ranging from 1 to max_value (C++)
Rcpp::cppFunction('NumericVector count_occurrences(NumericVector x, int max_value) {
  NumericVector result(max_value);
  for (int i = 0; i < x.length(); i++) { result[x(i)-1]++; }
  return result;
}')

#' Convert count statistics to rank (C++)
Rcpp::cppFunction('NumericVector count_to_rank(NumericVector x, int total_count) {
  NumericVector result(x.length());
  int cumulated_total = 0;
  for (int i = 0; i < x.length(); i++) {
    result[i] = (cumulated_total + (x(i)+1)*0.5) / total_count;
    cumulated_total += x[i];
  }
  return result;
}')

#' Convert bin number to rank (C++, in place operation)
Rcpp::cppFunction('void bin_to_rank(NumericVector bins, NumericVector rank_per_bin) {
  for (int i = 0; i < bins.length(); i++) { bins[i] = rank_per_bin[bins[i]-1]; }
}')

#' Compute neighbor voting from cell x cell correlation network
#'
#' voter_id is a cell x labels binary matrix indicating cell types
#' If left empty, cell types are assumed to be the column names of the network.
compute_votes_from_network <- function(network, voter_id = NULL) {
  if (is.null(voter_id)) {
    voter_id <- design_matrix(colnames(network))
  }
  return(network %*% voter_id / rowSums(network))
}

#' Compute AUROCs based on neighbor voting and candidate identities
#'
#' candidate_id is a binary matrix indicating the cell type of candidates
#' If left empty, cell types are assumed to be the row names of the vote matrix.
compute_aurocs <- function(votes, candidate_id = NULL) {
  if (is.null(candidate_id)) {
    positives <- design_matrix(rownames(votes))
  } else {
    positives <- candidate_id
  }
  n_positives <- colSums(positives)
  n_negatives <- nrow(positives) - n_positives
  sum_of_positive_ranks <- t(positives) %*% apply(abs(votes), MARGIN = 2, FUN = rank)
  result <- (sum_of_positive_ranks / n_positives - (n_positives+1)/2) / n_negatives
  return(result)
}

#' Transform a vector with cell_type labels into a binary matrix
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

#' Create a result matrix with one row and one column for each unique cell type
create_result_matrix <- function(cell_type) {
  unique_cell_type <- unique(cell_type)
  result <- matrix(0, nrow = length(unique_cell_type), ncol = length(unique_cell_type))
  rownames(result) <- unique_cell_type
  colnames(result) <- unique_cell_type
  return(result)
}

#' Return study id from a label in format 'study_id|cell_type'
get_study_id <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "\\|"), head, 1))
}

#' Return cell type from a label in format 'study_id|cell_type'
get_cell_type <- function(cluster_name) {
  return(sapply(strsplit(cluster_name, "\\|"), tail, 1))
}
