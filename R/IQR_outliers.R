#' This function computes interquartile range (IQR) criterion, which can be used
#' as a criterion for outlier detection
#'
#' @param DistVal numeric. vector of distribution of values
#' @param weightIRQ numeric. weighting factor applied to IRQ to define
#' lower and upper boundaries for outliers
#'
#' @return outlier_IQR numeric. band numbers of original sensor corresponding to S2
#' @importFrom stats IQR quantile
#' @export

IQR_outliers <- function(DistVal,weightIRQ = 1.5){
  range_IQR <- c(stats::quantile(DistVal, 0.25,na.rm=TRUE),
                 stats::quantile(DistVal, 0.75,na.rm=TRUE))
  iqr <- diff(range_IQR)
  outlier_IQR <- c(range_IQR[1]-weightIRQ*iqr,range_IQR[2]+weightIRQ*iqr)
  return(outlier_IQR)
}
