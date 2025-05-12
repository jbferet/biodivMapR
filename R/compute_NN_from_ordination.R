#' compute the nearest neighbors among kernels
#
#' @param mat_bc matrix of BC dissimilarity between the kernels excluded from Ordination (rows)
#' @param knn numeric. number of neighbors
#' @param pcoa_train numeric. PCoA coordinates computed from dissimilarity matrix
#
#' @return ordin_est estimated NMDS position based on nearest neighbors from NMDS
#' @export

compute_nn_from_ordination <- function(mat_bc, knn, pcoa_train) {
  # get nearest neighbors (coordinates and ID)
  nn <- apply(mat_bc, 1, FUN = sort,index.return = TRUE)
  ordin_est <- lapply(X = nn,
                      FUN = weighted_coords_nn,
                      knn = knn, pcoa_train = pcoa_train)
  ordin_est <- do.call(rbind,ordin_est)
  return(ordin_est)
}
