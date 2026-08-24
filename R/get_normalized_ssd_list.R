#' get spectral species distribution from a list for all clusters
#' including those with 0 abundance
#'
#' @param ssd_list dataframe. spectral species distribution
#' @param windows_id numeric. window id corresponding to ssd_list
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#
#' @return ssd_map normalized distribution for all clusters including those with 0 abundance
#' @importFrom stats kmeans
#' @export

get_normalized_ssd_list <- function(ssd_list, windows_id, nb_clusters, pcelim = 0.02){

  ssd_list <- lapply(X = ssd_list,FUN = table)
  ssd_mat <- do.call('rbind', lapply(X = ssd_list, FUN = get_normalized_ssd_mat,
                                     nb_clusters = nb_clusters, pcelim = pcelim))
  ssd_mat <- data.frame(ssd_mat)
  names(ssd_mat) <- seq_len(nb_clusters)
  ssd_mat$windows_id <- windows_id
  ssd_mat$iter <- names(ssd_list)
  return(ssd_mat)
}
