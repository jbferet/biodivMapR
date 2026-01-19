#' compute spectral species from inputdata
#
#' @param inputdata dataframe. information extracted from raster data
#' @param Kmeans_info list. obtained from function get_kmeans
#' @param selected_bands numeric. bands selected from inputdata
#' @param nbCPU numeric. Number of CPUs available
#
#' @return SSwindow results of kmeans per window
#' @importFrom rlist list.cbind
#' @export

get_spectralSpecies <- function(inputdata, Kmeans_info,
                                selected_bands = NULL, nbCPU = 1){
  if (is.null(selected_bands))
    selected_bands <- seq_len(ncol(inputdata))
  nb_iter <- length(Kmeans_info$Centroids)
  nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
  nb_pixels <- nrow(inputdata)
  # center reduce data
  inputdata_cr <- center_reduce(x = inputdata[selected_bands],
                                m = Kmeans_info$MinVal,
                                sig = Kmeans_info$Range)
  # inputdata_cr <- center_reduce(x = inputdata[,selected_bands],
  #                               m = Kmeans_info$MinVal,
  #                               sig = Kmeans_info$Range)
  # compute distance between each pixel and each centroid
  cluster_dist <- lapply(X = Kmeans_info$Centroids,
                         FUN = kmeans_iter,
                         inputvals = inputdata_cr)
  SSchunk <- rlist::list.cbind(cluster_dist)
  names(SSchunk) <- paste0('iter#',seq(1,ncol(SSchunk)))
  # # reshape distance into a matrix: all pixels from iteration 1, then all pixels from it2...
  # cluster_dist <- matrix(aperm(array(cluster_dist, c(nb_pixels, nb_clusters, nb_iter)), c(1, 3, 2)), nrow = nb_pixels * nb_iter)
  # # select closest cluster
  if (!is.null(inputdata$win_ID)) SSchunk$win_ID <- inputdata$win_ID
  return(SSchunk)
}
