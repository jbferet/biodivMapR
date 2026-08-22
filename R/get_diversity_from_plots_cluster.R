#' computes diversity metrics from validation plots
#'
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param rast_sample list. information extracted from raster
#' @param AttributeTable dataframe. corresponds to attributes of original vector file
#' @param Hill_order numeric. Hill order
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param selected_bands numeric. bands selected from input_rast
#' @param alpha_metrics character.
#' @param beta_metrics character.
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#'
#' @return list of spectral alpha and beta diversity metrics
#' @export

get_diversity_from_plots_cluster <- function(Kmeans_info, rast_sample,
                                             AttributeTable, Hill_order = 1,
                                             Beta_info = NULL,
                                             selected_bands = NULL,
                                             alpha_metrics = c('richness', 'shannon', 'simpson', 'hill'),
                                             beta_metrics = c('bray'),
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
  Attributes <- list()
  # Attributes <- data.frame(matrix(NA, ncol = ncol(AttributeTable),
  #                                 nrow = nbPlots_init))
  # names(Attributes) <- names(Attributes0)
  # Attributes[selPlots,] <- Attributes0

  windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
  windows_per_plot$win_ID <- list(SSValid$win_ID)

  Beta_info <- list()
  Beta_info$ssd <- list()
  for (j in seq_len(length(windows_per_plot$SSwindow_perCPU[[1]][[1]]))){
    ssdtmp <- list()
    for (i in seq_len(length(windows_per_plot$SSwindow_perCPU[[1]])))
      ssdtmp[[i]] <-  windows_per_plot$SSwindow_perCPU[[1]][[i]][[j]]
    ssdtmp <- lapply(X = ssdtmp,FUN = table)
    ssdtrain <- lapply(X = ssdtmp, FUN = get_normalized_ssd,
                       nb_clusters = nb_clusters, pcelim = pcelim)
    Beta_info$ssd[[j]] <- do.call('rbind', ssdtrain)

    # windows_per_plot <- split_chunk(SSchunk = aa[[1]], nbCPU = 1)

    # Beta_info$ssd[[j]] <- get_beta_from_ssd(ss_dist = Beta_info$ssd[[j]],
    #                                    beta_metrics =  beta_metrics,
    #                                    nb_clusters = nb_clusters,
    #                                    pcelim = pcelim)

  }
  diss_output <- TRUE
  alphabetaIdx_CPU <- lapply(X = windows_per_plot$SSwindow_perCPU,
                             FUN = alphabeta_window_list,
                             nb_clusters = nb_clusters,
                             alpha_metrics = alpha_metrics,
                             beta_metrics = beta_metrics,
                             Hill_order = Hill_order,
                             Beta_info = Beta_info,
                             diss_output = diss_output,
                             pcelim = pcelim)

  alphabetaIdx <- unlist(alphabetaIdx_CPU, recursive = FALSE)
  rm(alphabetaIdx_CPU)
  gc()
  # 7- reshape alpha diversity metrics
  IDwindow <- unlist(windows_per_plot$IDwindow_perCPU)
  res_shape_chunk <- list()

  for (idx in alpha_metrics){
    res_shape_chunk[[idx]] <- list()
    for (crit in c('mean', 'sd')){
      elem <- paste0(idx,'_',crit)
      restmp <- unlist(lapply(alphabetaIdx,'[[',elem))
      res_shape_chunk_tmp <- rep(x = NA,nbPlots_init)
      res_shape_chunk_tmp[IDwindow] <- restmp
      res_shape_chunk[[idx]][[crit]] <- res_shape_chunk_tmp
    }
  }


  mat_diss <- list()
  for (beta in beta_metrics)
    mat_diss[[beta]] <- matrix(data = NA, nrow = nbPlots_init, ncol = nbPlots_init)
  if ((!is.null(Beta_info) | diss_output) & !is.null(IDwindow)) {
    for (beta in beta_metrics){
      mat_diss0 <- do.call(rbind,lapply(alphabetaIdx,'[[',paste0('matrix_', beta)))
      mat_diss[[beta]][IDwindow,IDwindow] <- mat_diss0
      Attributes[[beta]] <- mat_diss[[beta]]
    }
  }

  Attributes$richness_mean <- res_shape_chunk$richness$mean
  Attributes$richness_sd <- res_shape_chunk$richness$sd
  Attributes$shannon_mean <- res_shape_chunk$shannon$mean
  Attributes$shannon_sd <- res_shape_chunk$shannon$sd
  Attributes$simpson_mean <- res_shape_chunk$simpson$mean
  Attributes$simpson_sd <- res_shape_chunk$simpson$sd
  Attributes$hill_mean <- res_shape_chunk$hill$mean
  Attributes$hill_sd <- res_shape_chunk$hill$sd
  return(list('specdiv' = Attributes))
}
