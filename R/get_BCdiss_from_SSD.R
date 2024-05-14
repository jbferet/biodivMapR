#' Computes BC dissimilarity for a list of spectral species distributions
#'
#' @param SSdist numeric. list of spectral species distribution
#' @param nbclusters numeric. number of clusters used in kmeans
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param p list. progressor object for progress bar
#'
#' @return list of spectral species distribution and corresponding  BC dissimilarity matrix corresponding to Mat1 and Mat2
#' @export

get_BCdiss_from_SSD <- function(SSdist, nbclusters, pcelim, p = NULL){
  SSD <- lapply(X = SSdist,FUN = table)
  SSD <- lapply(X = SSD,FUN = get_SSD_full,
                nbclusters = nbclusters, pcelim = pcelim)
  SSD <- do.call(rbind,SSD)
  SSDList <- list(SSD, SSD)
  MatBC <- compute_BCdiss(SSDList, pcelim)
  if (!is.null(p)) p()
  return(list('SSD'=SSD, 'MatBC' = MatBC))
}
