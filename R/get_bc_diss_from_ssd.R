#' Computes BC dissimilarity for a list of spectral species distributions
#'
#' @param ss_dist numeric. list of spectral species distribution
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param p list. progressor object for progress bar
#'
#' @return list of spectral species distribution and corresponding  BC
#' dissimilarity matrix corresponding to Mat1 and Mat2
#' @export

get_bc_diss_from_ssd <- function(ss_dist, nb_clusters, pcelim, p = NULL){
  ssd <- lapply(X = ss_dist,FUN = table)
  ssd <- lapply(X = ssd,FUN = get_normalized_ssd,
                nb_clusters = nb_clusters, pcelim = pcelim)
  ssd <- do.call(rbind,ssd)
  mat_bc <- dissUtils::diss(ssd, ssd, method = 'braycurtis')
  # ssd_list <- list(ssd, ssd)
  # mat_bc <- compute_bc_diss(ssd_list, pcelim)
  if (!is.null(p))
    p()
  return(list('SSD' = ssd, 'MatBC' = mat_bc))
}
