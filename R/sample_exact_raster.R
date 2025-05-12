#' sample exact number of pixels from a raster
#'
#' @param extent_area extent
#' @param nb_samples numeric. number of samples to be extracted
#' @param input_rast SpatRaster. raster to extract data from
#' @param input_mask SpatRaster. mask corresponding to raster to extract data from
#'
#' @return rast_sample dataframe. pixel/plot info extracted from input_rast
#' @importFrom sf st_sample st_as_sf
#' @importFrom terra vect buffer extract res
#' @importFrom stats runif
#' @export

sample_exact_raster <- function(extent_area = extent_area,
                                nb_samples = nb_samples,
                                input_rast = input_rast,
                                input_mask = input_mask){

  rast_sample <- sample_from_raster(extent_area = extent_area,
                                    nb_samples = nb_samples,
                                    input_rast = input_rast,
                                    input_mask = input_mask)
  return(rast_sample)
}
