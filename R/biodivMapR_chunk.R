#' apply biodivMapR (computes clusters + diversity metrics) to an image chunk
#'
#' @param blk list. rows and number of rows to read from
#' @param r_in list. path of rasters to read from
#' @param window_size numeric. window size for square plots
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#' @param fd_metrics character. list of functional metrics
#' @param selected_bands numeric. bands selected from input data
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#' @param p list. progressor object for progress bar
#'
#' @return Shannon index correspnding to the distribution
#' @import cli
#' @importFrom terra rast blocks readValues
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom dplyr filter %>%
#' @export

biodivMapR_chunk <- function(blk, r_in, window_size, Kmeans_info,
                             Beta_info = NULL, alpha_metrics = 'shannon',
                             Hill_order = 1, fd_metrics = NULL,
                             selected_bands = NULL, pcelim = 0.02, nbCPU = 1,
                             min_sun = 0.25, p = NULL){
  list_allidx <- c('richness', 'shannon', 'simpson', 'fisher', 'hill')
  # list_funct_idx <- fd_metrics <- c('FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')
  list_funct_idx <- c('FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')
  funct_idx_cpu <- NULL
  richness <- shannon <- simpson <- fisher <- hill <- list('mean' = NA,
                                                           'sd' = NA)
  # 1- read input files
  input_data <- res_shape_chunk <- list()
  nameVars <- c()
  for (fid in names(r_in)){
    input_data[[fid]] <- terra::readValues(r_in[[fid]], row = blk$row,
                                           nrows = blk$nrows, dataframe = TRUE)
    if (fid == 'mask')
      names(input_data[[fid]]) <- 'mask'
    if (dim(r_in[[fid]])[3]==1)
      nameVars <- c(nameVars, fid)
    if (dim(r_in[[fid]])[3]>1)
      nameVars <- names(input_data[[fid]])
  }
  if (is.null(selected_bands)){
    if ('mask'%in%nameVars)
      selected_bands <- seq_len(length(nameVars[-which(nameVars=='mask')]))
    if (!'mask'%in%nameVars)
      selected_bands <- seq_len(length(nameVars))
  }
  inputdata <- do.call(cbind, input_data)
  names(inputdata) <- nameVars
  inputdata$win_ID <- produce_win_ID(inputdata = inputdata, blk = blk,
                                     window_size = window_size)
  nbWindows <- max(inputdata$win_ID)
  # 2a- eliminate masked pixels
  if ('mask' %in% names(inputdata)){
    inputdata <- inputdata %>% dplyr::filter(inputdata$mask > 0)
    inputdata$mask <- NULL
  }
  # 2b- eliminate NA and inf
  inputdata <- clean_NAsInf(inputdata)
  if (nrow(inputdata)>0){
    # 3- eliminate windows with less than required number of sunlit/valid pixels
    inputdata <- get_sunlitwindows(inputdata = inputdata,
                                   pixperplot = window_size**2,
                                   min_sun = min_sun)
    # 4- convert pixel data to spectral species
    if (nrow(inputdata)>0){
      SSchunk <- get_spectralSpecies(inputdata = inputdata,
                                     Kmeans_info = Kmeans_info,
                                     selected_bands = selected_bands)
      # 5- split data chunk by window and by nbCPU to ensure parallel computing
      SSwindows_per_CPU <- split_chunk(SSchunk, nbCPU)
      # 6- compute diversity metrics
      nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
      if (nbCPU>1) {
        cl <- parallel::makeCluster(nbCPU)
        plan("cluster", workers = cl)
        alphabetaIdx_CPU <- future.apply::future_lapply(X = SSwindows_per_CPU$SSwindow_perCPU,
                                                        FUN = alphabeta_window_list,
                                                        nb_clusters = nb_clusters,
                                                        alpha_metrics = alpha_metrics,
                                                        Beta_info = Beta_info,
                                                        Hill_order = Hill_order,
                                                        pcelim = pcelim,
                                                        future.seed = TRUE)
        if (!is.null(fd_metrics)){
          inputdata <- cbind(center_reduce(x = inputdata[selected_bands],
                                           m = Kmeans_info$MinVal,
                                           sig = Kmeans_info$Range),
                             'win_ID' = inputdata$win_ID)
          windows_per_CPU <- split_chunk(inputdata, nbCPU)
          funct_idx_cpu <- future.apply::future_lapply(X = windows_per_CPU$SSwindow_perCPU,
                                                       FUN = functional_window_list,
                                                       fd_metrics = fd_metrics,
                                                       future.seed = TRUE)
        }
        parallel::stopCluster(cl)
        plan(sequential)
      } else {
        alphabetaIdx_CPU <- lapply(X = SSwindows_per_CPU$SSwindow_perCPU,
                                   FUN = alphabeta_window_list,
                                   nb_clusters = nb_clusters,
                                   alpha_metrics = alpha_metrics,
                                   Beta_info = Beta_info,
                                   Hill_order = Hill_order,
                                   pcelim = pcelim)
        if (!is.null(fd_metrics)){
          funct_idx_cpu <- lapply(X = windows_per_CPU$SSwindow_perCPU,
                                  FUN = functional_window_list,
                                  fd_metrics = fd_metrics)
        }
      }
      alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = FALSE)
      if (!is.null(fd_metrics))
        FunctionalIdx <- unlist(funct_idx_cpu, recursive = FALSE)
      rm(alphabetaIdx_CPU)
      rm(funct_idx_cpu)
      gc()
      # 7- reshape alpha diversity metrics
      IDwindow <- unlist(SSwindows_per_CPU$IDwindow_perCPU)
      for (idx in list_allidx){
        res_shape_chunk[[idx]] <- list()
        for (crit in c('mean', 'sd')){
          elem <- paste0(idx,'_',crit)
          restmp <- unlist(lapply(alphabetaIdx,'[[',elem))
          res_shape_chunk_tmp <- rep(x = NA,nbWindows)
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
          res_shape_chunk_tmp <- rep(x = NA,nbWindows)
          res_shape_chunk_tmp[IDwindow] <- restmp
          res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                           nrow = ceiling(blk$nrows/window_size),
                                           byrow = TRUE)
        } else {
          res_shape_chunk_tmp <- rep(x = NA,nbWindows)
          res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                           nrow = ceiling(blk$nrows/window_size),
                                           byrow = TRUE)
        }
      }
    } else {
      IDwindow <- NULL
      for (idx in list_allidx){
        res_shape_chunk[[idx]] <- list()
        for (crit in c('mean', 'sd')){
          elem <- paste0(idx,'_',crit)
          res_shape_chunk_tmp <- rep(x = NA,nbWindows)
          res_shape_chunk[[idx]][[crit]] <- matrix(res_shape_chunk_tmp,
                                                   nrow = ceiling(blk$nrows/window_size),
                                                   byrow = TRUE)
        }
      }
      for (idx in list_funct_idx){
        res_shape_chunk[[idx]] <- list()
        res_shape_chunk_tmp <- rep(x = NA,nbWindows)
        res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                         nrow = ceiling(blk$nrows/window_size),
                                         byrow = TRUE)
      }
    }
  } else {
    IDwindow <- NULL
    for (idx in list_allidx){
      res_shape_chunk[[idx]] <- list()
      for (crit in c('mean', 'sd')){
        elem <- paste0(idx,'_',crit)
        res_shape_chunk_tmp <- rep(x = NA,nbWindows)
        res_shape_chunk[[idx]][[crit]] <- matrix(res_shape_chunk_tmp,
                                                 nrow = ceiling(blk$nrows/window_size),
                                                 byrow = TRUE)
      }
    }
    for (idx in list_funct_idx){
      res_shape_chunk[[idx]] <- list()
      res_shape_chunk_tmp <- rep(x = NA,nbWindows)
      res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
                                       nrow = ceiling(blk$nrows/window_size),
                                       byrow = TRUE)
    }
  }
  # 8- reshape beta diversity metrics
  dimPCO <- 3
  if (!is.null(Beta_info))
    dimPCO <- ncol(Beta_info$BetaPCO$points)
  pcoa_bc <- matrix(data = NA, nrow = nbWindows, ncol = dimPCO)
  if (!is.null(Beta_info) & !is.null(IDwindow)) {
    pcoa_bc0 <- do.call(rbind,lapply(alphabetaIdx,'[[','PCoA_BC'))
    pcoa_bc[IDwindow,] <- pcoa_bc0
  }
  nb_rows <- ceiling(blk$nrows/window_size)
  nb_cols <- ceiling(nrow(pcoa_bc)/nb_rows)
  pcoa_bc <- aperm(array(data = c(pcoa_bc),dim = c(nb_cols,nb_rows,3)),c(2,1,3))
  if (!is.null(p))
    p()
  return(list('richness' = res_shape_chunk$richness,
              'shannon' = res_shape_chunk$shannon,
              'simpson' = res_shape_chunk$simpson,
              'fisher' = res_shape_chunk$fisher,
              'hill' = res_shape_chunk$hill,
              'FRic' = res_shape_chunk$FRic,
              'FEve' = res_shape_chunk$FEve,
              'FDiv' = res_shape_chunk$FDiv,
              'FDis' = res_shape_chunk$FDis,
              'FRaoq' = res_shape_chunk$FRaoq,
              'PCoA_BC' = pcoa_bc))
}
