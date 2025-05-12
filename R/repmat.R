#' R equivalent of repmat (matlab)
#'
#' @param X initial matrix
#' @param m nb of replications in row dimension
#' @param n nb of replications in column dimension
#'
#' @return matrix with values replicated and tiles
#' @export

repmat <- function(X, m, n) {
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  replicmat <- matrix(t(matrix(X, mx, nx * n)),
                      mx * m, nx * n, byrow = TRUE)
  return(replicmat)
}
