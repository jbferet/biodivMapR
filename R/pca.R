#' Function to perform PCA on a matrix
#'
#' @param X matrix to apply PCA on
#' @param type PCA (no rescale) or SPCA (rescale)
#'
#' @return list of PCA parameters (PCs from X, mean, eigenvectors and values)
#' @importFrom stats prcomp
#' @export

pca <- function(X, type = 'SPCA') {
  if (type == 'SPCA') scale = TRUE
  if (type == 'PCA') scale = FALSE
  modPCA <- stats::prcomp(X, scale = scale)
  return(modPCA)
}
