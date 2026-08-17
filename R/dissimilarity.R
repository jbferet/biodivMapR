# File: R/dissimilarity.R
# Description: R wrapper functions for dissimilarity metrics implemented in C++
# Author: Jean-Baptiste Feret
# Date: 2026-07-27

#' @title Compute Jaccard Dissimilarity
#' @description Computes the Jaccard dissimilarity between two matrices.
#' @param A A numeric matrix of size N x p (N observations, p species).
#' @param B A numeric matrix of size M x p (M observations, p species).
#' @return A dissimilarity matrix of size N x M.
#' @examples
#' A <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2, ncol = 3)
#' B <- matrix(c(0, 1, 1, 1, 0, 1), nrow = 2, ncol = 3)
#' jaccard_dissimilarity(A, B)
#' @useDynLib biodivMapR
#' @importFrom Rcpp evalCpp
#' @export
jaccard_dissimilarity <- function(A, B) {
  if (ncol(A) != ncol(B)) {
    stop("Matrices A and B must have the same number of columns (species).")
  }
  compute_jaccard_dissimilarity(A, B)
}

#' @title Compute Jaccard Turnover
#' @description Computes the Jaccard turnover between two matrices.
#' @param A A numeric matrix of size N x p.
#' @param B A numeric matrix of size M x p.
#' @return A turnover matrix of size N x M.
#' @examples
#' A <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2, ncol = 3)
#' B <- matrix(c(0, 1, 1, 1, 0, 1), nrow = 2, ncol = 3)
#' jaccard_turnover(A, B)
#' @useDynLib biodivMapR
#' @importFrom Rcpp evalCpp
#' @export
jaccard_turnover <- function(A, B) {
  if (ncol(A) != ncol(B)) {
    stop("Matrices A and B must have the same number of columns (species).")
  }
  compute_jaccard_turnover(A, B)
}

#' @title Compute Sorensen Dissimilarity
#' @description Computes the Sorensen dissimilarity between two matrices.
#' @param A A numeric matrix of size N x p.
#' @param B A numeric matrix of size M x p.
#' @return A dissimilarity matrix of size N x M.
#' @examples
#' A <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2, ncol = 3)
#' B <- matrix(c(0, 1, 1, 1, 0, 1), nrow = 2, ncol = 3)
#' sorensen_dissimilarity(A, B)
#' @useDynLib biodivMapR
#' @importFrom Rcpp evalCpp
#' @export
sorensen_dissimilarity <- function(A, B) {
  if (ncol(A) != ncol(B)) {
    stop("Matrices A and B must have the same number of columns (species).")
  }
  compute_sorensen_dissimilarity(A, B)
}

#' @title Compute Simpson Dissimilarity
#' @description Computes the Simpson dissimilarity between two matrices.
#' @param A A numeric matrix of size N x p.
#' @param B A numeric matrix of size M x p.
#' @return A dissimilarity matrix of size N x M.
#' @examples
#' A <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2, ncol = 3)
#' B <- matrix(c(0, 1, 1, 1, 0, 1), nrow = 2, ncol = 3)
#' simpson_dissimilarity(A, B)
#' @useDynLib biodivMapR
#' @importFrom Rcpp evalCpp
#' @export
simpson_dissimilarity <- function(A, B) {
  if (ncol(A) != ncol(B)) {
    stop("Matrices A and B must have the same number of columns (species).")
  }
  compute_simpson_dissimilarity(A, B)
}

#' @title Compute Bray-Curtis Dissimilarity
#' @description Computes the Bray-Curtis dissimilarity between two matrices.
#' @param A A numeric matrix of size N x p.
#' @param B A numeric matrix of size M x p.
#' @return A dissimilarity matrix of size N x M.
#' @examples
#' A <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2, ncol = 3)
#' B <- matrix(c(0, 1, 1, 1, 0, 1), nrow = 2, ncol = 3)
#' bray_curtis_dissimilarity(A, B)
#' @useDynLib biodivMapR
#' @importFrom Rcpp evalCpp
#' @export
bray_curtis_dissimilarity <- function(A, B) {
  if (ncol(A) != ncol(B)) {
    stop("Matrices A and B must have the same number of columns (species).")
  }
  compute_bray_curtis_dissimilarity(A, B)
}

#' @title Compute Bray-Curtis Turnover
#' @description Computes the Bray-Curtis turnover between two matrices.
#' @param A A numeric matrix of size N x p.
#' @param B A numeric matrix of size M x p.
#' @return A turnover matrix of size N x M.
#' @examples
#' A <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 2, ncol = 3)
#' B <- matrix(c(0, 1, 1, 1, 0, 1), nrow = 2, ncol = 3)
#' bray_curtis_turnover(A, B)
#' @useDynLib biodivMapR
#' @importFrom Rcpp evalCpp
#' @export
bray_curtis_turnover <- function(A, B) {
  if (ncol(A) != ncol(B)) {
    stop("Matrices A and B must have the same number of columns (species).")
  }
  compute_bray_curtis_turnover(A, B)
}
