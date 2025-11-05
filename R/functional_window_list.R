#' apply functional_window to a list of lists
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param fd_metrics list. functional diversity metrics: FRic, FEve, FDiv
#'
#' @return list of alpha and beta diversity metrics
#' @export

functional_window_list <- function(SSwindow,
                                   fd_metrics = c('FRic', 'FEve', 'FDiv')){
  functionalIdx <- lapply(X = SSwindow,
                          FUN = functional_window,
                          fd_metrics = fd_metrics)
  return(functionalIdx)
}
