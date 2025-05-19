#' center and reduce data matrix based on known mean and SD
#'
#' @param x numeric. data matrix (each column is centered/reduced)
#' @param m numeric. mean of each variable in the data matrix
#' @param sig numeric. SD of each variable in the data matrix
#'
#' @return x numeric. Centered matrix
#' @export

center_reduce <- function(x, m, sig) {
  for (i in seq_len(ncol(x)))
    x[, i] <- (x[, i] - m[i]) / sig[i]
  return(x)
}
