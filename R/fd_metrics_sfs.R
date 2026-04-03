#' compute functional metrics during sfs
#'
#' @param fd_metrics character
#' @param obs2optimize numeric. number of clusters used in kmeans
#' @param rast_val numeric
#' @param Kmeans_info list.
#' @param IDplot character
#' @param corrMethod character
#'
#' @return assesssed value and correlation
#' @export

fd_metrics_sfs <- function(fd_metrics, obs2optimize, rast_val, Kmeans_info, IDplot, corrMethod){
  if (!is.null(fd_metrics)){
    # center reduce data
    inputdata_cr <- center_reduce(x = rast_val,
                                  m = Kmeans_info$MinVal,
                                  sig = Kmeans_info$Range)
    inputdata_cr$ID <- IDplot
    inputdata_cr <- inputdata_cr %>% split(.$ID)
    inputdata_cr <- lapply(inputdata_cr,
                           function(x) data.frame(x[ , !(names(x) %in% 'ID')]))
    FunctDiv <- lapply(X = inputdata_cr,
                       FUN = get_functional_diversity,
                       fd_metrics = fd_metrics)
    FunctDiv <- data.frame('FRic' = unlist(lapply(FunctDiv, '[[',1)),
                           'FEve' = unlist(lapply(FunctDiv, '[[',2)),
                           'FDiv' = unlist(lapply(FunctDiv, '[[',3)))
    for (idx in fd_metrics){
      Assess[[idx]] <- FunctDiv[[idx]]
      CorrVal[[idx]]  <- cor.test(obs2optimize[[idx]],
                                  Assess[[idx]], method = corrMethod)$estimate
    }
  } else {
    Assess <- CorrVal <- list()
  }
  return(list('Assess' = Assess, 'CorrVal' = CorrVal))
}
