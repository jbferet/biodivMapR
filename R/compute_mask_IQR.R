#' compute mask based on interquartile range criterion applied on input rasters
#
#' @param input_raster_path character. path for input rasters.
#' @param output_mask_path character. path for output mask
#' @param input_mask_path character. path for optional input mask
#' @param weightIRQ numeric. weight to define SD range for IQR
#' @param filetype character. GDAL driver
#'
#' @return estimated NMDS position based on nearest neighbors from NMDS
#' @importFrom terra rast values names varnames writeRaster
#' @export

compute_mask_iqr <- function(input_raster_path,
                             output_mask_path,
                             input_mask_path = NULL,
                             weightIRQ = 3,
                             filetype = 'GTiff'){
  # initialize output mask: either input_mask, or only 1
  output_mask <- 1+0*terra::rast(input_raster_path[[1]])
  if (!is.null(input_mask_path)) {
    input_mask <- terra::rast(input_mask_path)
    output_mask <- input_mask
  } else if (is.null(input_mask_path)) {
    input_mask <- output_mask
    terra::values(output_mask)[which(is.na(terra::values(input_mask)))] <- 0
    names(input_mask) <- terra::varnames(input_mask) <- 'mask'
  }
  # update mask based on IQR computed for each raster
  for (si in input_raster_path){
    input_rast <- terra::rast(si)
    masked_raster <- terra::mask(x = input_rast,
                                 mask = input_mask)
    IQR <- IQR_SpatRaster(masked_raster, weightIRQ = 3)
    elim <- which(terra::values(masked_raster) < IQR$lowBound |
                    terra::values(masked_raster) > IQR$upBound)
    if (length(elim)>0)
      terra::values(output_mask)[elim] <- 0
  }
  # write updated mask compiling all IQR
  terra::writeRaster(x = output_mask, filename = output_mask_path,
                     overwrite = TRUE, filetype = filetype)
  return(invisible())
}
