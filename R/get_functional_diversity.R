#' get functional diversity metrics from dataframe
#' This function was inspired from FD package
#' @param spectraits numeric. dataframe containing species in rows and trait values in columns
#' @param fd_metrics character. Functional diversity metric
#' @param p list. progressor object for progress bar
#
#' @return fd_metricss
#' @importFrom fundiversity fd_fric fd_fdiv fd_feve fd_fdis fd_raoq
#' @importFrom future plan
#' @export

get_functional_diversity <- function(spectraits,
                                     fd_metrics = c('FRic', 'FEve', 'FDiv'),
                                     p = NULL){
  FRic <- FDiv <- FEve <- FDis <- FRaoq <- NA
  future::plan(future::sequential)
  if ('FRic' %in% fd_metrics){
    FRic <- try({
      fundiversity::fd_fric(spectraits)$FRic
    }, silent = TRUE)
    if ( "try-error" %in% class(FRic)){
      FRic <- NA
    }
  }

  if ('FDiv' %in% fd_metrics){
    FDiv <- try({
       fundiversity::fd_fdiv(spectraits)$FDiv
    }, silent = TRUE)
    if ( "try-error" %in% class(FDiv)){
      FDiv <- NA
    }
  }

  if ('FEve' %in% fd_metrics){
    FEve <- try({
      fundiversity::fd_feve(spectraits)$FEve
    }, silent = TRUE)
    if ( "try-error" %in% class(FEve)){
      FEve <- NA
    }
  }

  if ('FDis' %in% fd_metrics){
    FDis <- try({
      fundiversity::fd_fdis(spectraits)$FDis
    }, silent = TRUE)
    if ( "try-error" %in% class(FDis)){
      FDis <- NA
    }
  }

  if ('FRaoq' %in% fd_metrics){
    FRaoq <- try({
      fundiversity::fd_raoq(spectraits)$Q
    }, silent = TRUE)
    if ( "try-error" %in% class(FRaoq)){
      FRaoq <- NA
    }
  }

  if (!is.null(p))
    p()
  return(list('FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv,
              'FDis' = FDis, 'FRaoq' = FRaoq))
}
