#' computes kmeans for an iteration in biodivMapR
#'
#' @param Kmeans_Centroids numeric. coordinates of cluster centroids
#' @param inputvals numeric. chunk of input values to apply kmeans on
#'
#' @return SSchunk spectral species corresponding to the input chunk
#' @importFrom fields rdist
#' @export

kmeans_iter <- function(Kmeans_Centroids, inputvals){
  cluster_dist <- fields::rdist(inputvals, Kmeans_Centroids)
  SSchunk <- as.data.frame(max.col(-cluster_dist))
  names(SSchunk) <- 'cluster'
  rm(cluster_dist)
  return(SSchunk)
}
