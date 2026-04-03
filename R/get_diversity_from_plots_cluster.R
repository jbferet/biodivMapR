#' computes diversity metrics from validation plots
#'
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param rast_sample list. information extracted from raster
#' @param AttributeTable dataframe. corresponds to attributes of original vector file
#' @param Hill_order numeric. Hill order
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param selected_bands numeric. bands selected from input_rast
#' @param alpha_metrics character.
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#'
#' @return SpatVector including diversity metrics and BC dissimilarity for the plots
#' @importFrom dplyr group_split
#' @importFrom stats as.dist
#' @export

get_diversity_from_plots_cluster <- function(Kmeans_info, rast_sample,
                                             AttributeTable, Hill_order = 1,
                                             Beta_info = NULL,
                                             selected_bands = NULL,
                                             alpha_metrics = c('richness', 'shannon', 'simpson', 'hill'),
                                             pcelim = 0.02){

  nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
  win_ID <- NULL
  # get spectral species
  ssvect <- spectralspecies_per_rast_sample(Kmeans_info = Kmeans_info,
                                            rast_sample = rast_sample,
                                            selected_bands = selected_bands)
  SSValid <- ssvect$SSValid
  nbPlots_init <- length(unique(rast_sample$ID))
  selPlots <- seq_len(nbPlots_init)
  Attributes0 <- AttributeTable
  Attributes <- data.frame(matrix(NA, ncol = ncol(AttributeTable),
                                  nrow = nbPlots_init))
  names(Attributes) <- names(Attributes0)
  Attributes[selPlots,] <- Attributes0

  windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
  windows_per_plot$win_ID <- list(SSValid$win_ID)

  alphabetaIdx_CPU <- lapply(X = windows_per_plot$SSwindow_perCPU,
                             FUN = alphabeta_window_list,
                             nb_clusters = nb_clusters,
                             alpha_metrics = alpha_metrics,
                             Hill_order = Hill_order,
                             Beta_info = Beta_info,
                             pcelim = pcelim)

  alphabetaIdx <- unlist(alphabetaIdx_CPU, recursive = FALSE)
  rm(alphabetaIdx_CPU)
  gc()
  # 7- reshape alpha diversity metrics
  IDwindow <- unlist(windows_per_plot$IDwindow_perCPU)
  res_shapeChunk <- list()
  for (i in seq_len(10)) {
    restmp <- unlist(lapply(alphabetaIdx,'[[',i))
    res_shapeChunk[[i]] <- rep(x = NA,nbPlots_init)
    res_shapeChunk[[i]][IDwindow] <- restmp
  }

  Attributes$richness_mean <- res_shapeChunk[[1]]
  Attributes$richness_sd <- res_shapeChunk[[2]]
  Attributes$shannon_mean <- res_shapeChunk[[3]]
  Attributes$shannon_sd <- res_shapeChunk[[4]]
  Attributes$simpson_mean <- res_shapeChunk[[5]]
  Attributes$simpson_sd <- res_shapeChunk[[6]]
  Attributes$hill_mean <- res_shapeChunk[[9]]
  Attributes$hill_sd <- res_shapeChunk[[10]]
  return(list('specdiv' = Attributes))
}
