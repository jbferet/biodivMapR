#' Compute kmeans from random subset of pixels extracted from an image and a list
#' of values for k
#'
#' @param input_rast character. Path and name of the image to be processed.
#' @param input_mask character. Path and name of the mask corresponding to the image
#' @param SelectBands numeric. bands selected from input_rast
#' @param nbClust_list numeric. number of clusters used in kmeans
#' @param nbIter numeric. nb of iterations averaged to compute diversity indices
#' @param maxPixel_kmeans numeric. max number of pixels to extract for kmeans
#' @param algorithm character. algorithm used in the kmeans clustering
#' @param nbCPU numeric. Number of CPUs available
#' @param progressbar boolean. set true for progress bar during clustering
#'
#' @return Kmeans_info
#' @importFrom dplyr select all_of
#' @export

explore_kmeans <- function(input_rast,
                           input_mask = NULL,
                           SelectBands = NULL,
                           nbClust_list = 50,
                           nbIter = 20,
                           maxPixel_kmeans = 1e5,
                           algorithm = 'Hartigan-Wong',
                           nbCPU = 1, progressbar = F){

  # sample data from image
  Pix_Per_Iter <- define_pixels_per_iter(input_rast = input_rast,
                                         input_mask = input_mask,
                                         nbPix = maxPixel_kmeans,
                                         nbIter = nbIter)
  # define raster extent where to randomly sample square plots
  extent_area <- get_raster_extent(input_rast[[1]])
  # sample plots for initialization of beta diversity
  nbSamples <- Pix_Per_Iter*nbIter
  rast_sample <- sample_from_raster(extent_area = extent_area,
                                    nbSamples = nbSamples,
                                    input_rast = input_rast,
                                    input_mask = input_mask)
  if (is.null(SelectBands)) SelectBands <- seq_len(dim(rast_sample)[2])
  rast_sample <- rast_sample %>% select(all_of(SelectBands))
  # multi-thread
  registerDoFuture()
  cl <- parallel::makeCluster(nbCPU)
  plan("cluster", workers = cl)

  get_kmeans_list <- function() {
    foreach(nbclust = nbClust_list) %dopar% {
      Kmeans_info <- get_kmeans(rast_sample = rast_sample,
                                nbIter = nbIter,
                                nbclusters = nbclust,
                                algorithm = algorithm,
                                progressbar = F)
      return(Kmeans_info)
    }
  }
  Kmeans_info <- get_kmeans_list()
  plan(sequential)
  return(Kmeans_info)
}
