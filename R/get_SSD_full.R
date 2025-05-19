#' get spectral species distribution for all clusters
#
#' @param ssd dataframe. spectral species distribution
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#
#' @return SSwindow results of kmeans per window
#' @importFrom stats kmeans
#' @export

get_ssd_full <- function(ssd, nb_clusters, pcelim = 0.02){
  ssd_map <- 0*vector(length = nb_clusters)
  keep_ss <- which(ssd >= pcelim * sum(ssd))
  cluster_id <- as.numeric(names(ssd))[keep_ss]
  if (length(cluster_id)>0)
    ssd_map[cluster_id] <- ssd[keep_ss]
  # normalization
  ssd_map <- matrix(data = ssd_map/sum(ssd_map), nrow = 1)
  return(ssd_map)
}
