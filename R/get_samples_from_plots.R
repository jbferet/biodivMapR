#' extracts samples defined by a vector layer from a raster
#'
#' @param x spatVector vector layer
#' @param y spatRaster raster
#'
#' @return res
#' @importFrom terra vect extract
#' @export
#' 
get_samples_from_plots <- function(x, y){
  x <- terra::vect(x)
  res <- terra::extract(x = y, y = x, raw = TRUE, ID = FALSE)
  res <- c(unlist(res))
  return(res)
}
