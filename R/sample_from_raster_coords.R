#' sample set of pixels defined by row and col from raster data
#'
#' @param input_rast SpatRaster. raster to extract data from
#' @param xy data.frame coordinates.
#' @param input_mask SpatRaster. mask corresponding to raster to extract data from
#'
#' @return rast_sample dataframe. pixel/plot info extracted from input_rast
#' @importFrom sf st_sample st_as_sf
#' @importFrom terra vect buffer extract res
#' @importFrom stats runif
#' @export

sample_from_raster_coords <- function(input_rast, xy, input_mask = NULL){

  cell <- terra::cellFromRowCol(object = input_rast,
                                row = xy$row,
                                col = xy$col)
  # get coordinates for pixels to sample
  coords <- terra::xyFromCell(object = input_rast, cell = cell)
  rast_sample <- sample_raster(input_rast = input_rast, pix2extract = coords)
  # account for mask if provided
  if (!is.null(input_mask)){
    mask_sample <- terra::extract(x = input_mask, y = coords)
    sel <- which(mask_sample[[1]]==1)
    rast_sample <- rast_sample[sel,]
    xy <- xy[sel,]
  }
  rast_sample <- clean_NAsInf(rast_sample)
  return(list('DataSubset' = rast_sample, 'nbPix2Sample' = nrow(rast_sample),"coordPix"=xy))
}

