#' Function to perform MNF
#'
#' @param X numeric. matrix to apply MNF on
#' @param coordPix dataframe to compute noise, cf get_random_subset_from_image
#' @param retx boolean.
#'
#' @return results of MNF applied on matrix
#' @importFrom stats cov
#' @export
# used in noise

# TODO: faire 2 fonctions: mnf et mnf.subset, de même pour noise, noise.subset
mnf <- function(X, coordPix=NULL, retx=TRUE){
  if(any(is.na(X))){
    stop('Pixels with NAs found in X. Remove NA pixels before trying again.')
  }

  if(length(dim(X))>3)
    stop('X has more than 3 dimensions.')

  nz <- noise(X, coordPix)
  Xdim <- dim(X)

  if(is.null(coordPix) && length(dim(X))>2){
    X <- matrix(X[1:(Xdim[1]-1), 1:(Xdim[2]-1),], nrow = Xdim[1]*Xdim[2])
    nz <- matrix(nz, nrow = Xdim[1]*Xdim[2])
  }

  Xc = scale(X, center = T, scale = F)

  covNoise <- stats::cov(nz)
  covXc <- stats::cov(Xc)
  eig <- eigen(solve(covNoise)%*%covXc)
  colnames(eig$vectors) = paste0('PC', 1:ncol(eig$vectors))
  modMNF <- list(sdev = sqrt(eig$values), rotation = eig$vectors,
                 center = colMeans(X), scale = FALSE)
  attr(modMNF, 'class') <- 'prcomp'
  #   eig_pairs = tofsims:::EigenDecompose(covXc, covNoise, 1, nrow(covNoise))
  #   vord = order(Re(eig_pairs$eigval), decreasing = T)
  #   eig_pairs$eigval = Re(eig_pairs$eigval)[vord]
  #   eig_pairs$eigvec = Re(eig_pairs$eigvec[, vord])
  #   modMNF = list(rotation=eig_pairs$eigvec,
  #                 sdev=sqrt(eig_pairs$eigval),
  #                 center=colMeans(X),
  #                 scale=FALSE)
  if(retx==T)
    modMNF$x= array(Xc %*% modMNF$rotation, dim = Xdim)

  return(modMNF)
}
