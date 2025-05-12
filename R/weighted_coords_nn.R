#' Compute weighted coordinates of a spatial unit based on nearest neighbors
#' used during PCoA
#
#' @param nn list. coordinates & ID of nearest neighbors from pcoa_train matrix
#' @param knn number of neighbors
#' @param pcoa_train PCoA coordinates of reference samples
#
#' @return estimated NMDS position based on nearest neighbors from NMDS
#' @export

weighted_coords_nn <- function(nn, knn, pcoa_train) {

  # get distance and ID of nn samples
  dist_nn <- nn$x[seq_len(knn)]
  id_nn <- nn$ix[seq_len(knn)]
  # final location weighted by location of nn
  # if exact same location as nearest neighbor
  if (dist_nn[1]==0){
    mds_pos <- pcoa_train[id_nn[1], ]
  } else {
    # total dissimilarity from k nearest neighbors to weight
    dist_tot <- 1/sum(1/dist_nn)
    if (ncol(pcoa_train)>1){
      mds_pos <- colSums((dist_tot/dist_nn) * pcoa_train[id_nn, ])
    } else {
      mds_pos <- sum((dist_tot/dist_nn) * pcoa_train[id_nn, ])
    }
  }
  ordin_est <- mds_pos
  return(ordin_est)
}
