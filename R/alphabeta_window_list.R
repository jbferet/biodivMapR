#' apply alphabeta_window to a list of lists
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nbclusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#'
#' @return list of alpha and beta diversity metrics
#' @export

alphabeta_window_list <- function(SSwindow, nbclusters, Beta_info = NULL,
                                  alphametrics = 'shannon', pcelim = 0.02){
  alphabetaIdx <- lapply(X = SSwindow,
                         FUN = alphabeta_window,
                         nbclusters = nbclusters,
                         Beta_info = Beta_info,
                         alphametrics = alphametrics,
                         pcelim = pcelim)
  return(alphabetaIdx)
}
