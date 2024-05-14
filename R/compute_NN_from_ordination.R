#' compute the nearest neighbors among kernels
#
#' @param MatBC matrix of BC dissimilarity between the kernels excluded from Ordination (rows)
#' @param knn numeric. number of neighbors
#' @param PCoA_train numeric. PCoA coordinates computed from dissimilarity matrix
#
#' @return Ordin_est estimated NMDS position based on nearest neighbors from NMDS
#' @export

compute_NN_from_ordination <- function(MatBC, knn, PCoA_train) {
  # get nearest neighbors (coordinates and ID)
  NN <- apply(MatBC, 1, FUN = sort,index.return = TRUE)
  Ordin_est <- lapply(X = NN,
                      FUN = WeightedCoordsNN,
                      knn = knn, PCoA_train = PCoA_train)
  Ordin_est <- do.call(rbind,Ordin_est)
  return(Ordin_est)
}
