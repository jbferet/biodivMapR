#' Compute kmeans from random subset of pixels extracted from an image
#'
#' @param input_rast character. Path and name of the image to be processed.
#' @param output_dir character. Path for output directory
#' @param input_mask character. Path and name of the mask corresponding to the image
#' @param selected_bands numeric. bands selected from input_rast
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param nb_iter numeric. nb of iterations averaged to compute diversity indices
#' @param Kmeans_info_save character. path where to save Kmeans_info
#' @param Kmeans_info_read character. path where to read Kmeans_info
#' @param nb_samples_alpha numeric. max number of pixels to extract for kmeans
#' @param algorithm character. algorithm used in the kmeans clustering
#' @param nbCPU numeric. Number of CPUs available
#' @param verbose boolean. set true for messages
#' @param progressbar boolean. set true for progress bar during clustering
#'
#' @return Kmeans_info
#' @export

init_kmeans <- function(input_rast,
                        output_dir,
                        input_mask = NULL,
                        selected_bands = NULL,
                        nb_clusters = 50,
                        nb_iter = 10,
                        Kmeans_info_save = NULL,
                        Kmeans_info_read = NULL,
                        nb_samples_alpha = 1e5, algorithm = 'Hartigan-Wong',
                        nbCPU = 1, verbose = TRUE, progressbar = TRUE){

  # if Kmeans_info_read directs towards RData: read the variable if exists
  if (!is.null(Kmeans_info_read)){
    if (file.exists(Kmeans_info_read))
      load(Kmeans_info_read)
    if (!file.exists(Kmeans_info_read)){
      print_error_message('Kmeans_info_file_missing')
      Kmeans_info_read <- NULL
    }
  }
  # if Kmeans_info_read == NULL: compute kmeans
  if (is.null(Kmeans_info_read)){
    # 2- sample data from PCA image
    Pix_Per_Iter <- define_pixels_per_iter(input_rast = input_rast,
                                           input_mask = input_mask,
                                           nb_pix = nb_samples_alpha,
                                           nb_iter = nb_iter)
    # define raster extent where to randomly sample square plots
    extent_area <- get_raster_extent(input_rast[[1]])
    # sample plots for initialization of beta diversity
    nb_samples <- Pix_Per_Iter*nb_iter
    if (verbose)
      message('sampling pixels to compute spectral species')
    rast_sample <- sample_from_raster(extent_area = extent_area,
                                      nb_samples = nb_samples,
                                      input_rast = input_rast,
                                      input_mask = input_mask)

    Kmeans_info <- init_kmeans_samples(rast_sample = rast_sample,
                                       output_dir = output_dir,
                                       selected_bands = selected_bands,
                                       nb_clusters = nb_clusters,
                                       nb_iter = nb_iter,
                                       Kmeans_info_save = Kmeans_info_save,
                                       algorithm = algorithm, nbCPU = nbCPU,
                                       verbose = verbose,
                                       progressbar = progressbar)
  }
  return(Kmeans_info)
}
