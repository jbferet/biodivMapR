#' This function computes interquartile range (IQR) for a SpatRaster
#'
#' @param input_rast SpatRaster.
#' @param weightIRQ numeric. weighting factor applied to IRQ to define
#' lower and upper boundaries for outliers
#'
#' @return outlier_IQR numeric. band numbers of original sensor corresponding to S2
#' @importFrom terra global
#' @export

IQR_SpatRaster <- function(input_rast,weightIRQ = 3){
  IQ <- terra::global(x = input_rast, fun= quantile, na.rm = TRUE)
  range_IQR <- data.frame('quartile1' = IQ[['X25.']],
                          'quartile3' = IQ[['X75.']])
  iqr <- range_IQR$quartile3 - range_IQR$quartile1
  outlier_IQR <- data.frame('lowBound' = range_IQR$quartile1 - weightIRQ*iqr,
                            'upBound' = range_IQR$quartile3 + weightIRQ*iqr)
  return(outlier_IQR)
}
