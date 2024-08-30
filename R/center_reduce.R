#' center and reduce data matrix based on known mean and SD
#'
#' @param X numeric. data matrix (each column is centered/reduced)
#' @param m numeric. mean of each variable in the data matrix
#' @param sig numeric. SD of each variable in the data matrix
#'
#' @return X numeric. Centered matrix
#' @export

center_reduce <- function(X, m, sig) {
  for (i in seq_len(ncol(X))) X[, i] <- (X[, i] - m[i]) / sig[i]
  return(X)
}
