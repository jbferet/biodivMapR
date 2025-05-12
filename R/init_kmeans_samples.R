#' Compute kmeans from random subset of pixels extracted from an image
#'
#' @param rast_sample data frame containing samples to use
#' @param output_dir character. Path for output directory
#' @param selected_bands numeric. bands selected from input_rast
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param nb_iter numeric. nb of iterations averaged to compute diversity indices
#' @param Kmeans_info_save character. path where to save Kmeans_info
#' @param algorithm character. algorithm used in the kmeans clustering
#' @param nbCPU numeric. Number of CPUs available
#' @param verbose boolean. set true for messages
#' @param progressbar boolean. set true for progress bar during clustering
#'
#' @return Kmeans_info
#' @importFrom dplyr select all_of
#' @export

init_kmeans_samples <- function(rast_sample,
                                output_dir,
                                selected_bands = NULL,
                                nb_clusters = 50,
                                nb_iter = 10,
                                Kmeans_info_save = NULL,
                                algorithm = 'Hartigan-Wong',
                                nbCPU = 1, verbose = TRUE, progressbar = TRUE){

  rast_sample <- clean_NAsInf(rast_sample)
  if (is.null(selected_bands)) selected_bands <- seq_len(dim(rast_sample)[2])
  rast_sample <- rast_sample %>% select(all_of(selected_bands))
  # 3- PERFORM KMEANS FOR EACH ITERATION & DEFINE SPECTRAL SPECIES
  if (verbose)
    message("perform k-means clustering for each subset and define centroids")
  Kmeans_info <- get_kmeans(rast_sample = rast_sample,
                            nb_iter = nb_iter,
                            nb_clusters = nb_clusters, algorithm = algorithm,
                            nbCPU = nbCPU)
  if (is.null(Kmeans_info_save))
    Kmeans_info_save <- file.path(output_dir,'Kmeans_info.RData')
  save(Kmeans_info, file = Kmeans_info_save)
  return(Kmeans_info)
}
