#'save diversity maps as raster data once moving window process done
#'
#' @param rast_obj spatRaster.
#' @param window_size numeric. window size for square plots
#'
#' @return rast_crop
#' @importFrom terra ext res crop
#' @importFrom sf st_buffer
#' @export

crop_border_tile <- function(rast_obj, window_size){
  # get extent
  ext_rast <- terra::ext(x = rast_obj)
  # define bbox
  bbox <- preprocS2::bbox_to_poly(x = sf::st_bbox(ext_rast),
                                  crs = terra::crs(rast_obj))
  # create buffer
  bbox_buff <- sf::st_buffer(x = bbox,
                             dist = -(1+(window_size/2))*terra::res(rast_obj)[1])

  rast_crop <- terra::crop(x = rast_obj, y = bbox_buff)
  return(rast_crop)
}
