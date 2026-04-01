#' compute spectral species from a raster
#'
#' @param input_rast spatRaster. data raster
#' @param input_mask spatRaster. mask raster
#' @param selected_bands numeric. bands selected from input data
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param output_dir path where to save outputs
#' @param filetype character. gdal driver for output raster
#' @param overwrite boolean.
#'
#' @return none
#' @export

apply_spectral_species <- function(input_rast, input_mask = NULL,
                                   selected_bands = NULL, Kmeans_info,
                                   output_dir, filetype = 'GTiff',
                                   overwrite = TRUE){

  output_raster_full_name <- file.path(output_dir, 'spectral_species.tiff')
  if (FALSE %in% file.exists(output_raster_full_name) | overwrite){
    input_rast <- input_rast[[selected_bands]]
    # prepare to read input raster data
    r_in <- input_rast
    # if a mask file is provided
    selected_bands <- seq_len(dim(input_rast)[3])
    if (!is.null(input_mask))
      r_in[['mask']] <- input_mask
    inputdata <- as.data.frame(terra::values(r_in))
    if (!is.null(input_mask)){
      sel <- which(inputdata$mask==1)
    } else {
      sel <- seq_len(nrow(inputdata))
    }
    inputdata <- inputdata[sel,]
    ss_rast <- NA*input_rast[[1]]
    names(ss_rast) <- 'spectral_species'
    if (length(inputdata)>0){
      ss_tile <- get_spectralSpecies(inputdata = inputdata,
                                     Kmeans_info = Kmeans_info,
                                     selected_bands = selected_bands)
      ss_rast[sel] <- ss_tile[[1]]
    }
    terra::writeRaster(x = ss_rast, filename = output_raster_full_name,
                       filetype = filetype, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW"))
  }
  return(output_raster_full_name)
}

