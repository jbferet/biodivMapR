#' apply biodivMapR (computes clusters + diversity metrics) to an image chunk
#'
#' @param blk list. rows and number of rows to read from
#' @param r_in list. path of rasters to read from
#' @param window_size numeric. window size for square plots
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#' @param beta_metrics list. beta diversity metrics
#' @param fd_metrics character. list of functional metrics
#' @param selected_bands numeric. bands selected from input data
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#' @param p list. progressor object for progress bar
#'
#' @return alpha, beta and functional spectral diversity metrics corresponding to the chunk
#' @import cli
#' @importFrom terra rast blocks readValues
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom dplyr filter %>%
#' @export

biodivMapR_chunk <- function(
    blk, r_in, window_size, Kmeans_info, Beta_info = NULL,
    alpha_metrics = 'shannon', Hill_order = 1, beta_metrics = 'bray',
    fd_metrics = NULL, selected_bands = NULL, pcelim = 0.02, nbCPU = 1,
    min_sun = 0.25, p = NULL){

  # define list of alpha, beta and functional diversity metrics available
  list_alpha_idx <- c('richness', 'shannon', 'simpson', 'hill')
  list_beta_idx <- c('bray', 'brayturn', 'simpson_diss', 'jaccard', 'jaccardturn', 'sorensen')
  list_funct_idx <- c('FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')

  # initialize variables
  funct_idx_cpu <- NULL
  richness <- shannon <- simpson <- hill <- list('mean' = NA, 'sd' = NA)
  res_shape_chunk <- list()

  # 1- read input files
  l_in <- get_input_chunk(r_in, blk, selected_bands, window_size)
  input_data <- l_in$input_data
  selected_bands <- l_in$selected_bands
  nb_windows <- l_in$nb_windows
  l_in <- NULL

  if (nrow(input_data)>0){
    # 3- eliminate windows with less than required number of sunlit/valid pixels
    input_data <- get_sunlitwindows(inputdata = input_data,
                                   pixperplot = window_size**2,
                                   min_sun = min_sun)
    # 4- convert pixel data to spectral species
    if (nrow(input_data)>0){
      SSchunk <- get_spectralSpecies(inputdata = input_data,
                                     Kmeans_info = Kmeans_info,
                                     selected_bands = selected_bands)
      # 5- split data chunk by window and by nbCPU to ensure parallel computing
      SSwindows_per_CPU <- split_chunk(SSchunk, nbCPU)
      # 6- compute diversity metrics
      nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
      if (nbCPU>1) {
        cl <- parallel::makeCluster(nbCPU)
        parallel::clusterEvalQ(cl, {library(biodivMapR)})
        with(plan("cluster", workers = cl), local = TRUE)
        funct <- future.apply::future_lapply
        future.seed = TRUE
      } else {
        funct <- lapply
        future.seed = NULL
      }
      alphabetaIdx <- funct(X = SSwindows_per_CPU$SSwindow_perCPU,
                            FUN = alphabeta_chunk,
                            nb_clusters = nb_clusters,
                            alpha_metrics = alpha_metrics,
                            beta_metrics = beta_metrics,
                            Beta_info = Beta_info,
                            Hill_order = Hill_order,
                            pcelim = pcelim,
                            future.seed = future.seed,
                            future.chunk.size = NULL,
                            future.scheduling = 1)

      if (!is.null(fd_metrics)){
        input_data <- cbind(center_reduce(x = input_data[selected_bands],
                                         m = Kmeans_info$MinVal,
                                         sig = Kmeans_info$Range),
                           'win_ID' = input_data$win_ID)
        windows_per_CPU <- split_chunk(input_data, nbCPU)
        funct_idx_cpu <- funct(X = windows_per_CPU$SSwindow_perCPU,
                               FUN = functional_window_list,
                               fd_metrics = fd_metrics,
                               future.seed = future.seed)
      }
      if (nbCPU>1) {
        parallel::stopCluster(cl)
        plan(sequential)
      }

      if (!is.null(fd_metrics))
        FunctionalIdx <- unlist(funct_idx_cpu, recursive = FALSE)
      rm(funct_idx_cpu)
      gc()
      # 7- reshape alpha diversity metrics
      IDwindow <- unlist(SSwindows_per_CPU$IDwindow_perCPU)
      for (idx in list_alpha_idx){
        res_shape_chunk[[idx]] <- list()
        for (crit in c('mean', 'sd')){
          elem <- paste0(idx,'_',crit)
          restmp <- unlist(lapply(alphabetaIdx,'[[',elem))
          res_shape_chunk_tmp <- rep(x = NA,nb_windows)
          res_shape_chunk_tmp[IDwindow] <- restmp
          res_shape_chunk[[idx]][[crit]] <- matrix(res_shape_chunk_tmp,
                                                   nrow = ceiling(blk$nrows/window_size),
                                                   byrow = TRUE)
        }
      }
      for (idx in list_funct_idx){
        res_shape_chunk[[idx]] <- list()
        if (!is.null(fd_metrics)){
          restmp <- unlist(lapply(FunctionalIdx,'[[',idx))
          res_shape_chunk_tmp <- rep(x = NA,nb_windows)
          res_shape_chunk_tmp[IDwindow] <- restmp
          res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                           nrow = ceiling(blk$nrows/window_size),
                                           byrow = TRUE)
        } else {
          res_shape_chunk_tmp <- rep(x = NA,nb_windows)
          res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                           nrow = ceiling(blk$nrows/window_size),
                                           byrow = TRUE)
        }
      }
    } else {
      IDwindow <- NULL
      for (idx in list_alpha_idx){
        res_shape_chunk[[idx]] <- list()
        for (crit in c('mean', 'sd')){
          elem <- paste0(idx,'_',crit)
          res_shape_chunk_tmp <- rep(x = NA,nb_windows)
          res_shape_chunk[[idx]][[crit]] <- matrix(res_shape_chunk_tmp,
                                                   nrow = ceiling(blk$nrows/window_size),
                                                   byrow = TRUE)
        }
      }
      for (idx in list_funct_idx){
        res_shape_chunk[[idx]] <- list()
        res_shape_chunk_tmp <- rep(x = NA,nb_windows)
        res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                         nrow = ceiling(blk$nrows/window_size),
                                         byrow = TRUE)
      }
    }
  } else {
    IDwindow <- NULL
    for (idx in list_alpha_idx){
      res_shape_chunk[[idx]] <- list()
      for (crit in c('mean', 'sd')){
        elem <- paste0(idx,'_',crit)
        res_shape_chunk_tmp <- rep(x = NA,nb_windows)
        res_shape_chunk[[idx]][[crit]] <- matrix(res_shape_chunk_tmp,
                                                 nrow = ceiling(blk$nrows/window_size),
                                                 byrow = TRUE)
      }
    }
    for (idx in list_funct_idx){
      res_shape_chunk[[idx]] <- list()
      res_shape_chunk_tmp <- rep(x = NA,nb_windows)
      res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                       nrow = ceiling(blk$nrows/window_size),
                                       byrow = TRUE)
    }
  }
  # 8- reshape beta diversity metrics
  dimPCO <- 3
  if (!is.null(Beta_info))
    dimPCO <- ncol(Beta_info$beta_pco[[1]]$points)
  pcoa_diss <- list()
  for (beta in list_beta_idx)
    pcoa_diss[[beta]] <- matrix(data = NA, nrow = nb_windows, ncol = dimPCO)
  if (!is.null(Beta_info) & !is.null(IDwindow)) {
    for (beta in beta_metrics){
      pcoa_diss0 <- do.call(rbind,lapply(alphabetaIdx,'[[',paste0('pcoa_', beta)))
      pcoa_diss[[beta]][IDwindow,] <- pcoa_diss0
      nb_rows <- ceiling(blk$nrows/window_size)
      nb_cols <- ceiling(nrow(pcoa_diss[[beta]])/nb_rows)
      pcoa_diss[[beta]] <- aperm(array(data = c(pcoa_diss[[beta]]),dim = c(nb_cols,nb_rows,3)),c(2,1,3))
    }
  }
  if (!is.null(p))
    p()
  return(list('richness' = res_shape_chunk$richness,
              'shannon' = res_shape_chunk$shannon,
              'simpson' = res_shape_chunk$simpson,
              'hill' = res_shape_chunk$hill,
              'FRic' = res_shape_chunk$FRic,
              'FEve' = res_shape_chunk$FEve,
              'FDiv' = res_shape_chunk$FDiv,
              'FDis' = res_shape_chunk$FDis,
              'FRaoq' = res_shape_chunk$FRaoq,
              'pcoa_bray' = pcoa_diss[['bray']],
              'pcoa_brayturn' = pcoa_diss[['brayturn']],
              'pcoa_simpson_diss' = pcoa_diss[['simpson_diss']],
              'pcoa_jaccard' = pcoa_diss[['jaccard']],
              'pcoa_jaccardturn' = pcoa_diss[['jaccardturn']],
              'pcoa_sorensen' = pcoa_diss[['sorensen']]))
}
