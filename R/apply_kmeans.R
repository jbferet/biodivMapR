#' apply kmeans to information extracted from an image and corresponding to a window
#'
#' @param inputdata_window numeric. pixel data extracted from raster
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param SelectBands numeric. bands selected from input data
#'
#' @return Kmeans_info
#' @export

apply_kmeans <- function(inputdata_window, Kmeans_info, SelectBands = NULL){
  if (is.null(SelectBands)) SelectBands <- seq_len(dim(inputdata_window)[2])
  # get nb of partitions and nb of clusters based on Kmeans_info
  nbIter <- length(Kmeans_info$Centroids)
  nbclusters <- dim(Kmeans_info$Centroids[[1]])[1]
  nbPixels <- nrow(inputdata_window)
  # apply center reduction of raw raster data
  inputdata_cr <- center_reduce(X = inputdata_window[,SelectBands],
                                m = Kmeans_info$MinVal,
                                sig = Kmeans_info$Range)
  CentroidsArray <- do.call("rbind", Kmeans_info$Centroids)
  # compute distance between each pixel and each centroid
  cluster_dist <- rdist(inputdata_cr, CentroidsArray)
  # reshape distance into a matrix: all pixels from iteration 1, then all pixels from it2...
  cluster_dist <- matrix(aperm(array(cluster_dist, c(nbPixels, nbclusters, nbIter)),
                               c(1, 3, 2)), nrow = nbPixels * nbIter)
  ResDist <- as.data.frame(matrix(max.col(-cluster_dist), nrow = nbPixels))
  return(ResDist)
}
