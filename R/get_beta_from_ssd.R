#' Computes beta diversity for a list of spectral species distributions
#'
#' @param ss_dist numeric. list of spectral species distribution
#' @param beta_metrics character. name of beta dissimilarity metrics
#' Available options are "bray", "brayturn", "jaccard", "jaccardturn",
#' "sorensen", "simpson_diss"
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param p list. progressor object for progress bar
#' @return list of spectral species distribution and corresponding  BC
#' dissimilarity matrix corresponding to Mat1 and Mat2
#' @export

get_beta_from_ssd <- function(ss_dist, beta_metrics = 'bray', 
                              nb_clusters, pcelim, p = NULL){
  
  ssd <- lapply(X = ss_dist,FUN = table)
  ssd <- lapply(X = ssd,FUN = get_normalized_ssd,
                nb_clusters = nb_clusters, pcelim = pcelim)
  ssd <- do.call(rbind,ssd)
  mat_diss <- compute_dissimilarity(A = ssd, B = ssd, beta_metrics = beta_metrics)
  # mat_bc <- dissUtils::diss(ssd, ssd, method = 'braycurtis')
  # ssd_list <- list(ssd, ssd)
  # mat_bc <- compute_bc_diss(ssd_list, pcelim)
  if (!is.null(p))
    p()
  return(list('ssd' = ssd, 'mat_diss' = mat_diss))
}
