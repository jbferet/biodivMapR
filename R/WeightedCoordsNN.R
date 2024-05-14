#' Compute the weighted coordinates of a spatial unit based on nearest neighbors used during PCoA
#
#' @param NN list. coordinates and ID of nearest neighbors from PCoA_train matrix
#' @param knn number of neighbors
#' @param PCoA_train PCoA coordinates of reference samples
#
#' @return estimated NMDS position based on nearest neighbors from NMDS
#' @export

WeightedCoordsNN <- function(NN, knn, PCoA_train) {

  # get distance and ID of NN samples
  DistNN <- NN$x[1:knn]
  IdNN <- NN$ix[1:knn]
  # final location weighted by location of NN
  # if exact same location as nearest neighbor
  if (DistNN[1]==0){
    MDSpos <- PCoA_train[IdNN[1], ]
  } else {
    # total dissimilarity from k nearest neighbors to weight
    Dist_Tot <- 1/sum(1/DistNN)
    if (ncol(PCoA_train)>1){
      MDSpos <- colSums((Dist_Tot/DistNN) * PCoA_train[IdNN, ])
    } else {
      MDSpos <- sum((Dist_Tot/DistNN) * PCoA_train[IdNN, ])
    }
  }
  Ordin_est <- MDSpos
  return(Ordin_est)
}
