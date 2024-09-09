#' get functional diversity metrics from dataframe
#' This function was inspired from FD package
#' @param spectraits numeric. dataframe containing species in rows and trait values in columns
#' @param FDmetric character. Functional diversity metric
#' @param p list. progressor object for progress bar
#
#' @return FDmetrics
#' @importFrom fundiversity fd_fric fd_fdiv fd_feve fd_fdis fd_raoq
#' @export

get_functional_diversity <- function(spectraits,
                                     FDmetric = c('FRic', 'FEve', 'FDiv'),
                                     p = NULL){
  FRic <- FDiv <- FEve <- FDis <- FRaoq <- NA
  if ('FRic' %in% FDmetric) FRic <- fundiversity::fd_fric(spectraits)$FRic
  if ('FDiv' %in% FDmetric) FDiv <- fundiversity::fd_fdiv(spectraits)$FDiv
  if ('FEve' %in% FDmetric) FEve <- fundiversity::fd_feve(spectraits)$FEve
  if ('FDis' %in% FDmetric) FDis <- fundiversity::fd_fdis(spectraits)$FDis
  if ('FRaoq' %in% FDmetric) FRaoq <- fundiversity::fd_raoq(spectraits)$Q
  if (!is.null(p)){p()}
  return(list('FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv, 'FDis' = FDis, 'FRaoq' = FRaoq))
}
