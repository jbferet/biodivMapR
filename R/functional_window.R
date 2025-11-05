#' compute functional diversity metrics from pixel data corresponding to
#' features extracted from a window
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param fd_metrics list. alpha diversity metrics: richness, shannon, simpson
#' @param p list. progressor object for progress bar
#'
#' @return list of alpha and beta diversity metrics
#' @importFrom stats sd
#' @export

functional_window <- function(SSwindow, fd_metrics = c('FRic', 'FEve', 'FDiv'),
                              p = NULL){
  fmetrics <- data.frame('FRic'= NA, 'FEve'= NA, 'FDiv'= NA,
                         'FDis' = NA, 'FRaoq' = NA)
  if (nrow(SSwindow)>ncol(SSwindow))
    fmetrics <- get_functional_diversity(spectraits = SSwindow,
                                         fd_metrics = fd_metrics)
  if (!is.null(p))
    p()
  return(fmetrics)
}
