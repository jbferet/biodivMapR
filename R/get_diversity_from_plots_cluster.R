#' computes diversity metrics from validation plots
#'
#' @param input_rast SpatRaster
#' @param validation_vect SpatVector
#' @param Hill_order numeric. Hill order
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param input_mask SpatRaster
#' @param fd_metrics character.
#' @param selected_bands numeric. bands selected from input_rast
#' @param rast_sample dataframe.
#' @param AttributeTable dataframe.
#' @param alpha_metrics character.
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param getBeta boolean. set true if computation of beta required
#' @param verbose boolean. set true for messages
#' @param p object
#'
#' @return SpatVector including diversity metrics and BC dissimilarity for the plots
#' @importFrom dplyr group_split
#' @importFrom stats as.dist
#' @export

get_diversity_from_plots_cluster <- function(input_rast, validation_vect,
                                             Hill_order = 1,
                                             Kmeans_info, Beta_info = NULL,
                                             input_mask  = NULL, fd_metrics = NULL,
                                             selected_bands = NULL,
                                             rast_sample = NULL, AttributeTable = NULL,
                                             alpha_metrics = c('richness', 'shannon', 'simpson', 'hill'),
                                             min_sun = 0.25, pcelim = 0.02, nbCPU = 1,
                                             getBeta = TRUE, verbose = FALSE,
                                             p = NULL){
  win_ID <- NULL
  # get nb_iter and nb_clusters
  nb_iter <- length(Kmeans_info$Centroids)
  nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
  # read vector data
  if (inherits(validation_vect,
               what = 'SpatVectorCollection') & is.null(rast_sample)){
    SSValid <- Attributes <- list()
    nbPlots_init <- 0
    for (ind_vect in seq_len(length(validation_vect))){
      ssvect <- spectralspecies_per_polygon(SpatVector = validation_vect[[ind_vect]],
                                            input_rast = input_rast,
                                            fd_metrics = fd_metrics,
                                            selected_bands = selected_bands,
                                            input_mask = input_mask,
                                            Kmeans_info = Kmeans_info,
                                            rast_sample = rast_sample,
                                            AttributeTable = AttributeTable,
                                            min_sun = min_sun)
      if (!is.null(ssvect$SSValid)){
        SSValid[[ind_vect]] <- ssvect$SSValid
        Attributes[[ind_vect]] <- ssvect$AttributeTable
        SSValid[[ind_vect]]$win_ID <- SSValid[[ind_vect]]$win_ID + nbPlots_init
        Attributes[[ind_vect]]$ID_biodivMapR <- Attributes[[ind_vect]]$ID_biodivMapR + nbPlots_init
        nbPlots_init <- nbPlots_init + length(validation_vect[[ind_vect]])
      }
    }
    SSValid <- do.call(rbind,SSValid)
    Attributes <- do.call(rbind,Attributes)
  } else if (inherits(validation_vect, what = 'SpatVector') | (!is.null(rast_sample))){
    ssvect <- spectralspecies_per_polygon(SpatVector = validation_vect,
                                          input_rast = input_rast,
                                          input_mask = input_mask,
                                          fd_metrics = fd_metrics,
                                          Kmeans_info = Kmeans_info,
                                          selected_bands = selected_bands,
                                          rast_sample = rast_sample,
                                          AttributeTable = AttributeTable,
                                          min_sun = min_sun)
    SSValid <- ssvect$SSValid
    if (inherits(validation_vect, what = 'SpatVector')) {
      nbPlots_init <- length(validation_vect)
      nbPlots <- nrow(ssvect$AttributeTable)
      selPlots <- ssvect$AttributeTable$ID_biodivMapR
    } else if (!is.null(rast_sample)) {
      nbPlots_init <- nbPlots <- length(unique(rast_sample$ID))
      selPlots <- seq_len(nbPlots_init)
    }
    Attributes0 <- ssvect$AttributeTable
    Attributes <- data.frame(matrix(NA, ncol = ncol(ssvect$AttributeTable),
                                    # nrow = nrow(ssvect$AttributeTable)))
                                    nrow = nbPlots_init))
    names(Attributes) <- names(Attributes0)
    Attributes[selPlots,] <- Attributes0
  }

  windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
  windows_per_plot$win_ID <- list(SSValid$win_ID)

  alphabetaIdx_CPU <- lapply(X = windows_per_plot$SSwindow_perCPU,
                             FUN = alphabeta_window_list,
                             nb_clusters = nb_clusters,
                             alpha_metrics = alpha_metrics,
                             Hill_order = Hill_order,
                             Beta_info = Beta_info,
                             pcelim = pcelim)

  alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = FALSE)
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
  if (!is.null(p))
    p()
  return(list('specdiv' = Attributes))
}
