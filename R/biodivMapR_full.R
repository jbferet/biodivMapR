#' computes diversity metrics from raster
#'
#' @param input_raster_path character. path for the input rasters
#' @param output_dir character. path for the output files
#' @param window_size numeric. window size for square plots
#' @param maxRows numeric. max number of rows in each block
#' @param Kmeans_info_save character. path where to save Kmeans_info
#' @param Kmeans_info_read character. path where to read Kmeans_info
#' @param Beta_info_save character. path where to save Beta_info
#' @param Beta_info_read character. path where to read Beta_info
#' @param input_mask_path character. path for mask file
#' @param nbclusters numeric. number of clusters used in kmeans
#' @param nbSamples_beta numeric. number of samples to compute beta diversity from
#' @param SelectBands numeric. bands selected from input_rast
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param Hill_order numeric. Hill order
#' @param FDmetric character. list of functional metrics
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param nbIter numeric. nb of iterations averaged to compute diversity indices
#' @param MinSun numeric. minimum amount of sunlit pixels in the plots
#' @param maxPixel_kmeans numeric. max number of pixels extracted for kmeans
#' @param dimPCoA numeric. number of dimensions of PCoA
#' @param verbose boolean. set true for messages
#' @param progressbar boolean. set true for progress bar during clustering
#' @param filetype character. driver for output diversity raster data
#'
#' @return Kmeans_info and Beta_info
#' @export

biodivMapR_full <- function(input_raster_path, output_dir, window_size,
                            maxRows = NULL,
                            Kmeans_info_save = NULL, Kmeans_info_read = NULL,
                            Beta_info_save = NULL, Beta_info_read = NULL,
                            input_mask_path = NULL, nbclusters = 50,
                            nbSamples_beta = 1000, SelectBands = NULL,
                            alphametrics = 'shannon', Hill_order = 1, FDmetric = NULL,
                            pcelim = 0.02, nbCPU = 1, nbIter = 20, MinSun = 0.25,
                            maxPixel_kmeans = 1e5, dimPCoA = 3, verbose = T,
                            progressbar = T, filetype = 'GTiff'){

  # read input rasters
  if (inherits(x = input_raster_path, what = 'character')){
    input_rast <- terra::rast(input_raster_path)
  } else if (inherits(x = input_raster_path, what = 'list')){
    input_rast <- lapply(input_raster_path,terra::rast)
  }
  input_mask <- NULL
  if (!is.null(input_mask_path)) input_mask <- terra::rast(input_mask_path)
  if (is.null(Kmeans_info_save)) Kmeans_info_save <- file.path(output_dir,'Kmeans_info.RData')
  if (is.null(Beta_info_save)) Beta_info_save <- file.path(output_dir,'Beta_info.RData')
  # compute kmeans from random subset of image
  Kmeans_info <- init_kmeans(input_rast = input_rast,
                             output_dir = output_dir,
                             maxPixel_kmeans = maxPixel_kmeans,
                             input_mask = input_mask,
                             SelectBands = SelectBands,
                             nbclusters = nbclusters,
                             Kmeans_info_save = Kmeans_info_save,
                             Kmeans_info_read = Kmeans_info_read,
                             nbCPU = nbCPU)

  # compute beta diversity for training data
  Beta_info <- init_PCoA(input_rast = input_rast,
                         input_mask = input_mask,
                         window_size = window_size,
                         nbSamples = nbSamples_beta,
                         Kmeans_info = Kmeans_info, SelectBands = SelectBands,
                         nbCPU = nbCPU,
                         Beta_info_save = Beta_info_save,
                         Beta_info_read = Beta_info_read)

  # compute alpha and beta diversity from raster data
  # input_rasters <- list('main' = input_raster_path,
  #                       'mask' = input_mask_path)
  ab_div_metrics <- get_raster_diversity(input_raster_path = input_raster_path,
                                         input_mask_path = input_mask_path,
                                         Kmeans_info = Kmeans_info,
                                         Beta_info = Beta_info,
                                         SelectBands = SelectBands,
                                         window_size = window_size,
                                         alphametrics = alphametrics,
                                         Hill_order = Hill_order,
                                         FDmetric = FDmetric,
                                         pcelim = pcelim,
                                         maxRows = maxRows, nbCPU = nbCPU,
                                         MinSun = MinSun)

  # save diversity metrics as raster data
  save_diversity_maps(ab_div_metrics = ab_div_metrics,
                      alphametrics = alphametrics,
                      Hill_order = Hill_order,
                      FDmetric = FDmetric,
                      input_rast = input_rast,
                      output_dir = output_dir,
                      window_size = window_size,
                      filetype = filetype)
  return(list('Kmeans_info' = Kmeans_info,
              'Beta_info' = Beta_info))
}
