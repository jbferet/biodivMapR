#' Compute kmeans from random subset of pixels extracted from an image and a list
#' of values for k
#'
#' @param input_rast character. Path and name of the image to be processed.
#' @param input_mask character. Path and name of the mask corresponding to the image
#' @param selected_bands numeric. bands selected from input_rast
#' @param nbClust_list numeric. number of clusters used in kmeans
#' @param nb_iter numeric. nb of iterations averaged to compute diversity indices
#' @param nb_samples_alpha numeric. max number of pixels to extract for kmeans
#' @param algorithm character. algorithm used in the kmeans clustering
#' @param nbCPU numeric. Number of CPUs available
#' @param progressbar boolean. set true for progress bar during clustering
#'
#' @return kmeans_info
#' @importFrom dplyr select all_of
#' @importFrom doFuture registerDoFuture
#' @export

explore_kmeans <- function(input_rast,
                           input_mask = NULL,
                           selected_bands = NULL,
                           nbClust_list = 50,
                           nb_iter = 10,
                           nb_samples_alpha = 1e5,
                           algorithm = 'Hartigan-Wong',
                           nbCPU = 1, progressbar = FALSE){

  # sample data from image
  Pix_Per_Iter <- define_pixels_per_iter(input_rast = input_rast,
                                         input_mask = input_mask,
                                         nb_pix = nb_samples_alpha,
                                         nb_iter = nb_iter)
  # define raster extent where to randomly sample square plots
  extent_area <- get_raster_extent(input_rast[[1]])
  # sample plots for initialization of beta diversity
  nb_samples <- Pix_Per_Iter*nb_iter
  rast_sample <- sample_from_raster(extent_area = extent_area,
                                    nb_samples = nb_samples,
                                    input_rast = input_rast,
                                    input_mask = input_mask)
  if (is.null(selected_bands))
    selected_bands <- seq_len(dim(rast_sample)[2])
  rast_sample <- rast_sample %>% select(all_of(selected_bands))
  # multi-thread
  registerDoFuture()
  cl <- parallel::makeCluster(nbCPU)
  with(plan("cluster", workers = cl), local = TRUE)	
  nbclust <- NULL
  get_kmeans_list <- function() {
    foreach(nbclust = nbClust_list) %dopar% {
      kmeans_info <- get_kmeans(rast_sample = rast_sample,
                                nb_iter = nb_iter,
                                nb_clusters = nbclust,
                                algorithm = algorithm,
                                progressbar = FALSE)
      return(kmeans_info)
    }
  }
  kmeans_info <- get_kmeans_list()
  plan(sequential)
  return(kmeans_info)
}
