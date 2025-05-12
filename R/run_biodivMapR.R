#' apply biodivMapR to an individual set of rasters (one site)
#'
#' @param input_raster_path character.
#' @param input_mask_path directory where spectral indices are for each plot
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param output_dir path where to save outputs
#' @param output_raster_name raster file names
#' @param selected_bands numeric. bands selected from input data
#' @param window_size numeric. window size for square plots
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param Hill_order numeric. Hill order
#' @param FDmetric character. list of functional metrics
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#' @param maxRows numeric. maximum number or rows to process
#' @param nbCPU numeric. Number of CPUs available
#' @param min_sun numeric. minimum proportion of sunlit pixels
#' @param filetype character. gdal driver name
#' @param moving_window boolean. should process be moving window (much longer)
#'
#' @return none
#' @export

run_biodivMapR <- function(input_raster_path, input_mask_path = NULL,
                           Kmeans_info, Beta_info,
                           output_dir, output_raster_name,
                           selected_bands = NULL, window_size,
                           alphametrics = 'shannon',
                           Hill_order = 1, FDmetric = NULL, pcelim = 0.02,
                           maxRows = NULL, nbCPU = 1, min_sun = 0.25,
                           filetype = 'GTiff',
                           moving_window = FALSE){

  # read input rasters
  if (inherits(x = input_raster_path, what = 'character'))
    input_rast <- terra::rast(input_raster_path)
  if (inherits(x = input_raster_path, what = 'list'))
    input_rast <- lapply(input_raster_path,terra::rast)

  if (!moving_window){
    ab_div_metrics <- get_raster_diversity_tile(input_raster_path = input_raster_path,
                                                input_mask_path = input_mask_path,
                                                Kmeans_info = Kmeans_info,
                                                Beta_info = Beta_info,
                                                selected_bands = selected_bands,
                                                window_size = window_size,
                                                alphametrics = alphametrics,
                                                Hill_order = Hill_order,
                                                FDmetric = FDmetric,
                                                pcelim = pcelim,
                                                maxRows = maxRows,
                                                nbCPU = nbCPU,
                                                min_sun = min_sun)
    # save diversity metrics as raster data
    save_diversity_maps_tile(input_raster_path = input_raster_path,
                             ab_div_metrics = ab_div_metrics,
                             alphametrics = alphametrics,
                             Hill_order = Hill_order,
                             FDmetric = FDmetric,
                             input_rast = input_rast,
                             output_dir = output_dir,
                             output_raster_name = output_raster_name,
                             window_size = window_size,
                             filetype = filetype)
  }
  if (moving_window){
    ab_div_metrics <- get_raster_diversity_mw(input_raster_path = input_raster_path,
                                              input_mask_path = input_mask_path,
                                              Kmeans_info = Kmeans_info,
                                              Beta_info = Beta_info,
                                              selected_bands = selected_bands,
                                              window_size = window_size,
                                              alphametrics = alphametrics,
                                              Hill_order = Hill_order,
                                              FDmetric = FDmetric,
                                              pcelim = pcelim,
                                              maxRows = maxRows, nbCPU = nbCPU,
                                              min_sun = min_sun)

    save_diversity_maps_mw(input_raster_path = input_raster_path,
                           ab_div_metrics = ab_div_metrics,
                           alphametrics = alphametrics,
                           Hill_order = Hill_order,
                           FDmetric = FDmetric,
                           input_rast = input_rast,
                           output_dir = output_dir,
                           output_raster_name = output_raster_name,
                           window_size = window_size,
                           filetype = filetype)
  }
  return()
}
