#' get functional diversity metrics from dataframe
#' This function was inspired from FD package
#' @param spectraits numeric. dataframe containing species in rows and trait values in columns
#' @param fd_metrics character. Functional diversity metric
#' @param p list. progressor object for progress bar
#
#' @return fd_metricss
#' @importFrom fundiversity fd_fric fd_fdiv fd_feve fd_fdis fd_raoq
#' @export

get_functional_diversity <- function(spectraits,
                                     fd_metrics = c('FRic', 'FEve', 'FDiv'),
                                     p = NULL){
  FRic <- FDiv <- FEve <- FDis <- FRaoq <- NA
  if ('FRic' %in% fd_metrics) FRic <- fundiversity::fd_fric(spectraits)$FRic
  if ('FDiv' %in% fd_metrics) FDiv <- fundiversity::fd_fdiv(spectraits)$FDiv
  if ('FEve' %in% fd_metrics) FEve <- fundiversity::fd_feve(spectraits)$FEve
  if ('FDis' %in% fd_metrics) FDis <- fundiversity::fd_fdis(spectraits)$FDis
  if ('FRaoq' %in% fd_metrics) FRaoq <- fundiversity::fd_raoq(spectraits)$Q
  if (!is.null(p)){p()}
  return(list('FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv,
              'FDis' = FDis, 'FRaoq' = FRaoq))
}
