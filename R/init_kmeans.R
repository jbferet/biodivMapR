#' Compute kmeans from random subset of pixels extracted from an image
#'
#' @param input_rast character. Path and name of the image to be processed.
#' @param output_dir character. Path for output directory
#' @param input_mask character. Path and name of the mask corresponding to the image
#' @param SelectBands numeric. bands selected from input_rast
#' @param nbclusters numeric. number of clusters used in kmeans
#' @param nbIter numeric. nb of iterations averaged to compute diversity indices
#' @param Kmeans_info_save character. path where to save Kmeans_info
#' @param Kmeans_info_read character. path where to read Kmeans_info
#' @param maxPixel_kmeans numeric. max number of pixels to extract for kmeans
#' @param algorithm character. algorithm used in the kmeans clustering
#' @param nbCPU numeric. Number of CPUs available
#' @param verbose boolean. set true for messages
#' @param progressbar boolean. set true for progress bar during clustering
#'
#' @return Kmeans_info
#' @importFrom dplyr select all_of
#' @export

init_kmeans <- function(input_rast,
                        output_dir,
                        input_mask = NULL,
                        SelectBands = NULL,
                        nbclusters = 50,
                        nbIter = 20,
                        Kmeans_info_save = NULL,
                        Kmeans_info_read = NULL,
                        maxPixel_kmeans = 1e5, algorithm = 'Hartigan-Wong',
                        nbCPU = 1, verbose = T, progressbar = T){

  # if Kmeans_info_read directs towards RData: read the variable if exists
  if (!is.null(Kmeans_info_read)){
    if (file.exists(Kmeans_info_read)) load(Kmeans_info_read)
    if (!file.exists(Kmeans_info_read)){
      print_error_message('Kmeans_info_file_missing')
      Kmeans_info_read = NULL
    }
  }
  # if Kmeans_info_read == NULL: compute kmeans
  if (is.null(Kmeans_info_read)){
    # 2- sample data from PCA image
    Pix_Per_Iter <- define_pixels_per_iter(input_rast = input_rast,
                                           input_mask = input_mask,
                                           nbPix = maxPixel_kmeans,
                                           nbIter = nbIter)
    # define raster extent where to randomly sample square plots
    extent_area <- get_raster_extent(input_rast[[1]])
    # sample plots for initialization of beta diversity
    nbSamples <- Pix_Per_Iter*nbIter
    if (verbose ==T)
      message('sampling pixels to compute spectral species')
    rast_sample <- sample_from_raster(extent_area = extent_area,
                                      nbSamples = nbSamples,
                                      input_rast = input_rast,
                                      input_mask = input_mask)

    Kmeans_info <- init_kmeans_samples(rast_sample = rast_sample,
                                       output_dir = output_dir,
                                       SelectBands = SelectBands,
                                       nbclusters = nbclusters, nbIter = nbIter,
                                       Kmeans_info_save = Kmeans_info_save,
                                       algorithm = algorithm, nbCPU = nbCPU,
                                       verbose = verbose, progressbar = progressbar)

    # if (is.null(SelectBands)) SelectBands <- seq_len(dim(rast_sample)[2])
    # rast_sample <- rast_sample %>% select(all_of(SelectBands))
    # # 3- PERFORM KMEANS FOR EACH ITERATION & DEFINE SPECTRAL SPECIES
    # if (verbose ==T)
    #   message("perform k-means clustering for each subset and define centroids")
    # Kmeans_info <- get_kmeans(rast_sample = rast_sample, nbIter = nbIter,
    #                           nbclusters = nbclusters, algorithm = algorithm,
    #                           nbCPU = nbCPU)
    # if (is.null(Kmeans_info_save))
    #   means_info_save <- file.path(output_dir,'Kmeans_info.RData')
    # save(Kmeans_info, file = Kmeans_info_save)
  }
  return(Kmeans_info)
}
