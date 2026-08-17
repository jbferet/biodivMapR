#' computes diversity metrics from raster
#'
#' @param input_raster_path character. path for the input rasters
#' @param input_mask_path character. path for mask file
#' @param output_dir character. path for the output files
#' @param selected_bands numeric. bands selected from input_rast
#' @param Kmeans_info_save character. path where to save Kmeans_info
#' @param Kmeans_info_read character. path where to read Kmeans_info
#' @param options list. includes options
#' - alpha_metrics list. alpha diversity metrics: richness, shannon, simpson
#' - Hill_order numeric. Hill order
#' - beta_metrics boolean. set TRUE to compute beta diversity
#' - fd_metrics character. list of functional metrics
#' - nb_samples_alpha numeric. max number of pixels extracted for kmeans
#' - nb_samples_beta numeric. number of samples to compute beta diversity
#' - nb_clusters numeric. number of clusters used in kmeans
#' - nb_iter numeric. nb of iterations averaged to compute diversity indices
#' - pcelim numeric. minimum proportion of pixels to consider spectral species
#' - maxRows numeric. max number of rows processed once by each CPU
#' - moving_window boolean. should diversity be computed on moving window?
#' - min_sun numeric. minimum amount of sunlit pixels in the plots
#' - dimPCoA numeric. number of dimensions of PCoA
#' - progressbar boolean. set true for progress bar during clustering
#' - filetype character. driver for output diversity raster data
#'
#' @return path for diversity_maps, Kmeans_info and Beta_info
#' @export
#'
spectral_species_full <- function(input_raster_path, input_mask_path = NULL,
                                  output_dir, selected_bands = NULL,
                                  Kmeans_info_save = NULL,
                                  Kmeans_info_read = NULL, options = NULL){

  # define options
  options <- set_options_biodivMapR(fun = 'spectral_species_full', options = options)
  nb_samples_alpha <- options$nb_samples_alpha
  nb_clusters <- options$nb_clusters
  nb_iter <- 1
  maxRows <- options$maxRows

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # read and/or arrange input rasters with terra 
  if (inherits(x = input_raster_path, what = 'list')){
    if (inherits(x = input_raster_path[[1]], what = 'character')){
      input_rast <- lapply(input_raster_path,terra::rast)
      input_rast <- terra::rast(input_rast)
    } else if (inherits(x = input_raster_path[[1]], what = 'SpatRaster')){
      input_rast <- terra::rast(input_raster_path)
    }
  } else if (inherits(x = input_raster_path, what = 'character')){
    input_rast <- terra::rast(input_raster_path)
  } else if (inherits(x = input_raster_path, what = 'SpatRaster')){
    input_rast <- input_raster_path
  }
  # # read input rasters
  # if (inherits(x = input_raster_path, what = 'character')){
  #   input_rast <- terra::rast(input_raster_path)
  # } else if (inherits(x = input_raster_path, what = 'list')){
  #   input_rast <- lapply(input_raster_path,terra::rast)
  # }
  input_mask <- NULL
  if (!is.null(input_mask_path))
    input_mask <- terra::rast(input_mask_path)
  if (is.null(Kmeans_info_save))
    Kmeans_info_save <- file.path(output_dir,'Kmeans_info.RData')
  # compute kmeans from random subset of image
  options <- set_options_biodivMapR(fun = 'init_kmeans', options = options)
  Kmeans_info <- init_kmeans(input_rast = input_rast,
                             output_dir = output_dir,
                             input_mask = input_mask,
                             selected_bands = selected_bands,
                             nb_clusters = nb_clusters,
                             Kmeans_info_save = Kmeans_info_save,
                             Kmeans_info_read = Kmeans_info_read,
                             nb_samples_alpha = nb_samples_alpha,
                             nb_iter = nb_iter, options = options)

  # apply clustering on raster
  rast_spectral_species <- apply_spectral_species(input_rast = input_rast,
                                                  input_mask = input_mask,
                                                  selected_bands = selected_bands,
                                                  Kmeans_info = Kmeans_info,
                                                  output_dir = output_dir,
                                                  filetype = 'GTiff',
                                                  overwrite = TRUE)

  return(list('rast_spectral_species' = rast_spectral_species,
              'Kmeans_info' = Kmeans_info))
}
