#' apply kmeans to information extracted from an image and corresponding to a
#' window
#'
#' @param inputdata_window numeric. pixel data extracted from raster
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param select_bands numeric. bands selected from input data
#'
#' @return Kmeans_info
#' @export

apply_kmeans <- function(inputdata_window, Kmeans_info, select_bands = NULL){
  if (is.null(select_bands)) select_bands <- seq_len(dim(inputdata_window)[2])
  # get nb of partitions and nb of clusters based on Kmeans_info
  nb_iter <- length(Kmeans_info$Centroids)
  nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
  nb_pixels <- nrow(inputdata_window)
  # apply center reduction of raw raster data
  inputdata_cr <- center_reduce(x = inputdata_window[,select_bands],
                                m = Kmeans_info$MinVal,
                                sig = Kmeans_info$Range)
  centroids_array <- do.call("rbind", Kmeans_info$Centroids)
  # compute distance between each pixel and each centroid
  cluster_dist <- rdist(inputdata_cr, centroids_array)
  # reshape distance into a matrix: pixels from it 1, then pixels from it2
  cluster_dist <- matrix(aperm(array(cluster_dist,
                                     c(nb_pixels, nb_clusters, nb_iter)),
                               c(1, 3, 2)), nrow = nb_pixels * nb_iter)
  res_dist <- as.data.frame(matrix(max.col(-cluster_dist), nrow = nb_pixels))
  return(res_dist)
}
