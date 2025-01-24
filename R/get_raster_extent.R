#' gets raster extent
#'
#' @param input_rast SpatRaster.
#'
#' @return extent_area as SpatVector
#' @importFrom terra ext vect
#' @export

get_raster_extent <- function(input_rast){
  # define data where to randomly sample square cells within the image and mask
  proj <- terra::crs(input_rast)
  x <- as.vector(terra::ext(input_rast))
  xCorners <- c(x[1L], x[2L], x[2L], x[1L])
  yCorners <- c(x[3L], x[3L], x[4L], x[4L])
  corners <- cbind(xCorners, yCorners)
  extent_area <- terra::vect(corners, type='polygon', crs=proj)
  return(extent_area)
}

