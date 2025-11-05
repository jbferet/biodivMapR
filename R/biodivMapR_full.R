#' computes diversity metrics from raster
#'
#' @param input_raster_path character. path for the input rasters
#' @param input_mask_path character. path for mask file
#' @param output_dir character. path for the output files
#' @param window_size numeric. window size for square plots
#' @param selected_bands numeric. bands selected from input_rast
#' @param Kmeans_info_save character. path where to save Kmeans_info
#' @param Kmeans_info_read character. path where to read Kmeans_info
#' @param Beta_info_save character. path where to save Beta_info
#' @param Beta_info_read character. path where to read Beta_info
#' @param nbCPU numeric. Number of CPUs available
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
biodivMapR_full <- function(input_raster_path, input_mask_path = NULL,
                            output_dir, window_size,
                            selected_bands = NULL, Kmeans_info_save = NULL,
                            Kmeans_info_read = NULL,  Beta_info_save = NULL,
                            Beta_info_read = NULL,  nbCPU = 1, options = NULL){

  # define options
  options <- set_options_biodivMapR(fun = 'biodivMapR_full', options = options)
  alpha_metrics <- options$alpha_metrics
  Hill_order <- options$Hill_order
  beta_metrics <- options$beta_metrics
  fd_metrics <- options$fd_metrics
  nb_samples_alpha <- options$nb_samples_alpha
  nb_samples_beta <- options$nb_samples_beta
  nb_clusters <- options$nb_clusters
  nb_iter <- options$nb_iter
  pcelim <- options$pcelim
  maxRows <- options$maxRows
  moving_window <- options$moving_window
  min_sun <- options$min_sun
  dimPCoA <- options$dimPCoA
  progressbar <- options$progressbar
  filetype <- options$filetype

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # read input rasters
  if (inherits(x = input_raster_path, what = 'character')){
    input_rast <- terra::rast(input_raster_path)
  } else if (inherits(x = input_raster_path, what = 'list')){
    input_rast <- lapply(input_raster_path,terra::rast)
  }
  input_mask <- NULL
  if (!is.null(input_mask_path))
    input_mask <- terra::rast(input_mask_path)
  if (is.null(Kmeans_info_save))
    Kmeans_info_save <- file.path(output_dir,'Kmeans_info.RData')
  if (is.null(Beta_info_save) & beta_metrics)
    Beta_info_save <- file.path(output_dir,'Beta_info.RData')
  # compute kmeans from random subset of image
  Kmeans_info <- init_kmeans(input_rast = input_rast,
                             output_dir = output_dir,
                             input_mask = input_mask,
                             selected_bands = selected_bands,
                             nb_clusters = nb_clusters,
                             Kmeans_info_save = Kmeans_info_save,
                             Kmeans_info_read = Kmeans_info_read,
                             nb_samples_alpha = nb_samples_alpha,
                             nb_iter = nb_iter,
                             nbCPU = nbCPU)

  # compute beta diversity for training data
  Beta_info <- NULL
  if (beta_metrics)
    Beta_info <- init_PCoA(input_rast = input_rast,
                           input_mask = input_mask,
                           window_size = window_size,
                           nb_samples = nb_samples_beta,
                           Kmeans_info = Kmeans_info,
                           selected_bands = selected_bands,
                           nbCPU = nbCPU, min_sun = min_sun,
                           pcelim = pcelim, dimPCoA = dimPCoA,
                           Beta_info_save = Beta_info_save,
                           Beta_info_read = Beta_info_read)

  # compute alpha and beta diversity from raster data
  # input_rasters <- list('main' = input_raster_path,
  #                       'mask' = input_mask_path)
  options(fundiversity.memoise = FALSE)

  if (!moving_window){
    message('compute raster diversity using moving window')
    message('please set "moving_window = FALSE" if this takes too much time')
    ab_div_metrics <- get_raster_diversity(input_raster_path = input_raster_path,
                                           input_mask_path = input_mask_path,
                                           Kmeans_info = Kmeans_info,
                                           Beta_info = Beta_info,
                                           selected_bands = selected_bands,
                                           window_size = window_size,
                                           alpha_metrics = alpha_metrics,
                                           Hill_order = Hill_order,
                                           fd_metrics = fd_metrics,
                                           pcelim = pcelim,
                                           maxRows = maxRows, nbCPU = nbCPU,
                                           min_sun = min_sun)

    # save diversity metrics as raster data
    diversity_maps <- save_diversity_maps(ab_div_metrics = ab_div_metrics,
                                          alpha_metrics = alpha_metrics,
                                          Hill_order = Hill_order,
                                          beta_metrics = beta_metrics,
                                          fd_metrics = fd_metrics,
                                          input_rast = input_rast,
                                          output_dir = output_dir,
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
                                              alpha_metrics = alpha_metrics,
                                              Hill_order = Hill_order,
                                              fd_metrics = fd_metrics,
                                              pcelim = pcelim,
                                              maxRows = maxRows, nbCPU = nbCPU,
                                              min_sun = min_sun)

    betanames <- 'beta_mw'
    alphanames <- paste0(alpha_metrics,'_mw')
    functionalname <- NULL
    if (!is.null(fd_metrics))
      functionalname <- paste0(fd_metrics,'_mw')
    output_raster_name <- as.list(c(betanames, alphanames, functionalname))
    names(output_raster_name) <- c('beta', alpha_metrics, fd_metrics)
    diversity_maps <- save_diversity_maps_mw(input_raster_path = input_raster_path,
                                             ab_div_metrics = ab_div_metrics,
                                             alpha_metrics = alpha_metrics,
                                             Hill_order = Hill_order,
                                             beta_metrics = beta_metrics,
                                             fd_metrics = fd_metrics,
                                             input_rast = input_rast,
                                             output_dir = output_dir,
                                             output_raster_name = output_raster_name,
                                             window_size = window_size,
                                             filetype = filetype)
  }
  return(list('diversity_maps' = diversity_maps,
              'Kmeans_info' = Kmeans_info,
              'Beta_info' = Beta_info))
}
