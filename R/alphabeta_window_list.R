#' apply alphabeta_window to a list of lists
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#'
#' @return list of alpha and beta diversity metrics
#' @export

alphabeta_window_list <- function(SSwindow, nb_clusters, Beta_info = NULL,
                                  alpha_metrics = 'shannon', Hill_order = 1,
                                  pcelim = 0.02){
  alphabeta_idx <- lapply(X = SSwindow,
                          FUN = alphabeta_window,
                          nb_clusters = nb_clusters,
                          Beta_info = Beta_info,
                          alpha_metrics = alpha_metrics,
                          Hill_order = Hill_order,
                          pcelim = pcelim)
  return(alphabeta_idx)
}
