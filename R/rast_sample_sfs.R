#' extract informaton from raster to prepare for SFS
#'
#' @param input_raster SpatRaster or list of SpatRaster
#' @param input_mask SpatRaster corresponding to mask
#' @param nb_pix numeric.
#' @param nb_iter numeric
#'
#' @return list
#'
#' @export

rast_sample_sfs <- function(input_raster, input_mask, nb_pix, nb_iter){
  Pix_Per_Iter <- define_pixels_per_iter(input_rast = input_raster,
                                         input_mask = input_mask,
                                         nb_pix = nb_pix,
                                         nb_iter = nb_iter)
  extent_area <- get_raster_extent(input_raster[[1]])
  nb_samples <- Pix_Per_Iter*nb_iter
  rast_sample <- sample_from_raster(extent_area = extent_area,
                                    nb_samples = nb_samples,
                                    input_rast = input_raster,
                                    input_mask = input_mask)
  return(rast_sample)
}
