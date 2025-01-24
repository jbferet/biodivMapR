#' computes shannon index from a distribution
#' (faster than version implemented in vegan package)
#'
#' @param Distrib Distribution
#'
#' @return Shannon index corresponding to the distribution
#' @export

get_Shannon <- function(Distrib) {
  Distrib <- Distrib / sum(Distrib, na.rm = TRUE)
  Distrib <- Distrib[which(!Distrib == 0)]
  shannon <- -1 * sum(Distrib * log(Distrib), na.rm = TRUE)
  return(shannon)
}
