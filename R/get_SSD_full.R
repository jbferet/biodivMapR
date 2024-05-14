#' get spectral species distribution for all clusters, even those with null
#' abundance
#
#' @param SSD dataframe. spectral species distribution
#' @param nbclusters numeric. number of clusters used in kmeans
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#
#' @return SSwindow results of kmeans per window
#' @importFrom stats kmeans
#' @export

get_SSD_full <- function(SSD, nbclusters, pcelim = 0.02){
  SSDMap0 <- 0*vector(length = nbclusters)
  KeepSS <- which(SSD >= pcelim * sum(SSD))
  ClusterID <- as.numeric(names(SSD))[KeepSS]
  if (length(ClusterID)>0) SSDMap0[ClusterID] <- SSD[KeepSS]
  # normalization
  SSDMap0 <- matrix(data = SSDMap0/sum(SSDMap0),nrow = 1)
  return(SSDMap0)
}
