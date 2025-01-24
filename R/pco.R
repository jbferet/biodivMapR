#' compute PCoA using the cmdscale function
#' original source: package labdsv https://rdrr.io/cran/labdsv/src/R/pco.R
#
#' @param dis matrix of BC dissimilarity between the kernels excluded from Ordination (rows)
#' @param k numeric. number of dimensions for PCoA
#
#' @return list
#' @importFrom stats cmdscale
#' @export
#'
pco <- function(dis, k=3){
  tmp <- stats::cmdscale(dis, k = k, eig = TRUE)
  class(tmp) <- c("dsvord","pco")
  tmp$type <- "PCO"
  return(tmp)
}
