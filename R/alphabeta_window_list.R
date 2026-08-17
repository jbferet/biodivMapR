#' apply alphabeta_window to a list of lists
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#' @param beta_metrics list. beta diversity metrics
#' @param diss_output boolean set TRUE to get dissimilarity matrices as outputs
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#'
#' @return list of alpha and beta diversity metrics
#' @export

alphabeta_window_list <- function(
    SSwindow,
    nb_clusters,
    Beta_info = NULL,
    alpha_metrics = 'shannon',
    Hill_order = 1,
    beta_metrics = 'bray',
    diss_output = FALSE,
    pcelim = 0.02){

  alphabeta_idx <- lapply(X = SSwindow,
                          FUN = alphabeta_window,
                          nb_clusters = nb_clusters,
                          Beta_info = Beta_info,
                          alpha_metrics = alpha_metrics,
                          Hill_order = Hill_order,
                          beta_metrics = beta_metrics,
                          diss_output = diss_output,
                          pcelim = pcelim)
  return(alphabeta_idx)
}
