#' apply functional_window to a list of lists
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param FDmetric list. functional diversity metrics: FRic, FEve, FDiv
#'
#' @return list of alpha and beta diversity metrics
#' @export

functional_window_list <- function(SSwindow,
                                   FDmetric = c('FRic', 'FEve', 'FDiv')){
  functionalIdx <- lapply(X = SSwindow,
                          FUN = functional_window,
                          FDmetric = FDmetric)
  return(functionalIdx)
}
