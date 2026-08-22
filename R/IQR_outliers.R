#' This function computes interquartile range (IQR) criterion, which can be used
#' as a criterion for outlier detection
#' outlier_range_iqr <- c(range_IQR$1]- weightIQR x iqr, range_IQR(2)+ weightIQR x iqr)
#' with range_IQR(1) the first quartile and range_IQR(2) the third quartile
#'
#' @param DistVal numeric. vector of distribution of values
#' @param weightIQR numeric. weighting factor applied to IQR to define
#' lower and upper boundaries for outliers
#'
#' @return outlier_range_iqr numeric. minimum and maximum values beyond which samples are considered outliers
#' @importFrom stats IQR quantile
#' @export

iqr_outliers <- function(DistVal,weightIQR = 1.5){
  range_IQR <- c(stats::quantile(DistVal, 0.25,na.rm=TRUE),
                 stats::quantile(DistVal, 0.75,na.rm=TRUE))
  iqr <- diff(range_IQR)
  outlier_range_iqr <- c(range_IQR[1]-weightIQR*iqr,range_IQR[2]+weightIQR*iqr)
  return(outlier_range_iqr)
}
