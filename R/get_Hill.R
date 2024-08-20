#' computes hill number from a distribution
#' The Hill numbers quantify biodiversity. The importance of the abundance
#' distribution increases with increasing Hill order.
#' For q=0, the Hill number is the richness, for q=1, it is the exponential
#' Shannon entropy and for q=2, it is the inverse Simpson index.
#' Note that the Hill order can also be a fraction, e.g. 0.5.
#' @param Distrib Distribution
#' @param q numeric. Hill order
#'
#' @return hill number corresponding to the distribution
#' @export

get_Hill <- function(Distrib, q = 1){
  Distrib <- Distrib / sum(Distrib, na.rm = TRUE)
  Distrib <- Distrib[which(!Distrib == 0)]
  if(q==0) {
    hill <- length(Distrib)
  } else if (q==1) {
    hill <- exp(get_Shannon(Distrib))
  } else {
    hill <- (sum(Distrib**q))**(1/(1-q))
  }
  return(hill)
}
