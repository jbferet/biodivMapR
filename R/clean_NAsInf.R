#' cleans dataframe from NAs and Inf values
#
#' @param df dataframe. df to be cleaned
#
#' @return df dataframe.
#' @importFrom stats na.omit
#' @export

clean_NAsInf <- function(df){
  # eliminate NAs
  df <- stats::na.omit(object = df)
  # eliminate infs
  Infval <- unlist(lapply(data.frame(lapply(df, is.infinite)), which))
  if (length(Infval)>0) df <- df[-Infval,]
  return(df)
}
