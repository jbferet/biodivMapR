#' applies results of ordination to full image based on nearest neighbors
#
#' @param x numeric.
#' @param centers numeric.
#' @param iter.max numeric.
#' @param nstart numeric.
#' @param algorithm character. algorithm used in the kmeans clustering
#' @param p function.
#
#' @return results of kmeans
#' @importFrom stats kmeans
#' @export

kmeans_progressr <- function(x, centers, iter.max, nstart,
                             algorithm = "Hartigan-Wong", p = NULL){
  res <- kmeans(x = x, centers = centers, iter.max = iter.max, nstart = nstart,
                algorithm = algorithm)
  if (!is.null(p)) p()
  return(res)
}
