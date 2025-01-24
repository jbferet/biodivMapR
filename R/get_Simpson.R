#' computes Simpson index from a distribution
#' (faster than version implemented in vegan package)
#'
#' @param Distrib Distribution
#'
#' @return Simpson index correspnding to the distribution
#' @export

get_Simpson <- function(Distrib) {
  Distrib <- Distrib / sum(Distrib, na.rm = TRUE)
  Simpson <- 1 - sum(Distrib * Distrib, na.rm = TRUE)
  return(Simpson)
}
