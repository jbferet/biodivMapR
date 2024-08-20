#' apply biodivMapR (computes clusters + diversity metrics) to an image chunk
#'
#' @param blk list. rows and number of rows to read from
#' @param r_in list. path of rasters to read from
#' @param window_size numeric. window size for square plots
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param Hill_order numeric. Hill order
#' @param FDmetric character. list of functional metrics
#' @param SelectBands numeric. bands selected from input data
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param MinSun numeric. minimum amount of sunlit pixels in the plots
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

biodivMapR_chunk <- function(blk, r_in, window_size, Kmeans_info, Beta_info = NULL,
                             alphametrics = 'shannon', Hill_order = 1, FDmetric = NULL,
                             SelectBands = NULL, pcelim = 0.02, nbCPU = 1,
                             MinSun = 0.25, p = NULL){
  list_allidx <- c('richness', 'shannon', 'simpson', 'fisher', 'hill')
  list_Funcidx <- c('FRic', 'FEve', 'FDiv')
  FunctionalIdx_CPU <- NULL
  richness <- shannon <- simpson <- fisher <- hill <- list('mean' = NA, 'sd' = NA)
  # 1- read input files
  input_data <- res_shapeChunk <- list()
  nameVars <- c()
  for (fid in names(r_in)){
    input_data[[fid]] <- terra::readValues(r_in[[fid]], row = blk$row,
                                           nrows = blk$nrows, dataframe = TRUE)
    if (fid == 'mask') names(input_data[[fid]]) <- 'mask'
    if (dim(r_in[[fid]])[3]==1) nameVars <- c(nameVars, fid)
    if (dim(r_in[[fid]])[3]>1) nameVars <- names(input_data[[fid]])
  }
  if (is.null(SelectBands)){
    if ('mask'%in%nameVars) SelectBands <- seq_len(length(nameVars[-which(nameVars=='mask')]))
    if (!'mask'%in%nameVars) SelectBands <- seq_len(length(nameVars))
  }
  inputdata <- do.call(cbind, input_data)
  names(inputdata) <- nameVars
  inputdata$win_ID <- produce_win_ID(inputdata = inputdata, blk = blk,
                                     window_size = window_size)
  nbWindows <- max(inputdata$win_ID)
  # 2a- eliminate masked pixels
  if ('mask' %in% names(inputdata)){
    inputdata <- inputdata %>% dplyr::filter(inputdata$mask == 1)
    inputdata$mask <- NULL
  }
  # 2b- eliminate NA and inf
  inputdata <- clean_NAsInf(inputdata)
  if (nrow(inputdata)>0){
    # 3- eliminate windows with less than required number of sunlit / valid pixels
    inputdata <- get_sunlitwindows(inputdata = inputdata,
                                   pixperplot = window_size**2,
                                   MinSun = MinSun)
    # 4- convert pixel data to spectral species
    if (nrow(inputdata)>0){
      SSchunk <- get_spectralSpecies(inputdata = inputdata,
                                     Kmeans_info = Kmeans_info,
                                     SelectBands = SelectBands)
      # 5- split data chunk by window and by nbCPU to ensure parallel computing
      whichID <- which(names(inputdata) =='win_ID')
      windows_per_CPU <- split_chunk(inputdata[c(SelectBands, whichID)], nbCPU)
      SSwindows_per_CPU <- split_chunk(SSchunk, nbCPU)
      # 6- compute diversity metrics
      nbclusters <- dim(Kmeans_info$Centroids[[1]])[1]
      if (nbCPU>1) {
        # plan(multisession, workers = nbCPU)
        cl <- parallel::makeCluster(nbCPU)
        plan("cluster", workers = cl)  ## same as plan(multisession, workers = nbCPU)
        alphabetaIdx_CPU <- future.apply::future_lapply(X = SSwindows_per_CPU$SSwindow_perCPU,
                                                        FUN = alphabeta_window_list,
                                                        nbclusters = nbclusters,
                                                        alphametrics = alphametrics,
                                                        Beta_info = Beta_info,
                                                        Hill_order = Hill_order,
                                                        pcelim = pcelim)
        if (!is.null(FDmetric)){
          FunctionalIdx_CPU <- future.apply::future_lapply(X = windows_per_CPU$SSwindow_perCPU,
                                                           FUN = functional_window_list,
                                                           FDmetric = FDmetric)
        }
        parallel::stopCluster(cl)
        plan(sequential)
      } else {
        alphabetaIdx_CPU <- lapply(X = SSwindows_per_CPU$SSwindow_perCPU,
                                   FUN = alphabeta_window_list,
                                   nbclusters = nbclusters,
                                   alphametrics = alphametrics,
                                   Beta_info = Beta_info,
                                   Hill_order = Hill_order,
                                   pcelim = pcelim)
        if (!is.null(FDmetric)){
          FunctionalIdx_CPU <- lapply(X = windows_per_CPU$SSwindow_perCPU,
                                      FUN = functional_window_list,
                                      FDmetric = FDmetric)
        }
      }
      alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = F)
      if (!is.null(FDmetric)) FunctionalIdx <- unlist(FunctionalIdx_CPU,recursive = F)
      rm(alphabetaIdx_CPU)
      rm(FunctionalIdx_CPU)
      gc()
      # 7- reshape alpha diversity metrics
      IDwindow <- unlist(SSwindows_per_CPU$IDwindow_perCPU)
      for (idx in list_allidx){
        res_shapeChunk[[idx]] <- list()
        for (crit in c('mean', 'sd')){
          elem <- paste0(idx,'_',crit)
          restmp <- unlist(lapply(alphabetaIdx,'[[',elem))
          res_shapeChunk_tmp <- rep(x = NA,nbWindows)
          res_shapeChunk_tmp[IDwindow] <- restmp
          res_shapeChunk[[idx]][[crit]] <- matrix(res_shapeChunk_tmp,
                                                  nrow = ceiling(blk$nrows/window_size),
                                                  byrow = T)
        }
      }
      for (idx in list_Funcidx){
        res_shapeChunk[[idx]] <- list()
        if (!is.null(FDmetric)){
          restmp <- unlist(lapply(FunctionalIdx,'[[',idx))
          res_shapeChunk_tmp <- rep(x = NA,nbWindows)
          res_shapeChunk_tmp[IDwindow] <- restmp
          res_shapeChunk[[idx]] <- matrix(res_shapeChunk_tmp,
                                          nrow = ceiling(blk$nrows/window_size),
                                          byrow = T)
        } else {
          res_shapeChunk_tmp <- rep(x = NA,nbWindows)
          res_shapeChunk[[idx]] <- matrix(res_shapeChunk_tmp,
                                          nrow = ceiling(blk$nrows/window_size),
                                          byrow = T)
        }
      }
    } else {
      IDwindow <- NULL
      for (idx in list_allidx){
        res_shapeChunk[[idx]] <- list()
        for (crit in c('mean', 'sd')){
          elem <- paste0(idx,'_',crit)
          res_shapeChunk_tmp <- rep(x = NA,nbWindows)
          res_shapeChunk[[idx]][[crit]] <- matrix(res_shapeChunk_tmp,
                                                  nrow = ceiling(blk$nrows/window_size),
                                                  byrow = T)
        }
      }
      for (idx in list_Funcidx){
        res_shapeChunk[[idx]] <- list()
        res_shapeChunk_tmp <- rep(x = NA,nbWindows)
        res_shapeChunk[[idx]] <- matrix(res_shapeChunk_tmp,
                                        nrow = ceiling(blk$nrows/window_size),
                                        byrow = T)
      }
    }
  } else {
    IDwindow <- NULL
    for (idx in list_allidx){
      res_shapeChunk[[idx]] <- list()
      for (crit in c('mean', 'sd')){
        elem <- paste0(idx,'_',crit)
        res_shapeChunk_tmp <- rep(x = NA,nbWindows)
        res_shapeChunk[[idx]][[crit]] <- matrix(res_shapeChunk_tmp,
                                                nrow = ceiling(blk$nrows/window_size),
                                                byrow = T)
      }
    }
    for (idx in list_Funcidx){
      res_shapeChunk[[idx]] <- list()
      res_shapeChunk_tmp <- rep(x = NA,nbWindows)
      res_shapeChunk[[idx]] <- matrix(res_shapeChunk_tmp,
                                      nrow = ceiling(blk$nrows/window_size),
                                      byrow = T)
    }
  }
  # 8- reshape beta diversity metrics
  dimPCO <- 3
  if (!is.null(Beta_info)) dimPCO <- ncol(Beta_info$BetaPCO$points)
  PCoA_BC <- matrix(data = NA, nrow = nbWindows, ncol = dimPCO)
  if (!is.null(Beta_info) & !is.null(IDwindow)) {
    PCoA_BC0 <- do.call(rbind,lapply(alphabetaIdx,'[[','PCoA_BC'))
    PCoA_BC[IDwindow,] <- PCoA_BC0
  }
  nbRows <- ceiling(blk$nrows/window_size)
  nbCols <- ceiling(nrow(PCoA_BC)/nbRows)
  PCoA_BC <- aperm(array(data = c(PCoA_BC),dim = c(nbCols,nbRows,3)),c(2,1,3))
  if (!is.null(p)) p()
  return(list('richness' = res_shapeChunk$richness,
              'shannon' = res_shapeChunk$shannon,
              'simpson' = res_shapeChunk$simpson,
              'fisher' = res_shapeChunk$fisher,
              'hill' = res_shapeChunk$hill,
              'FRic' = res_shapeChunk$FRic,
              'FEve' = res_shapeChunk$FEve,
              'FDiv' = res_shapeChunk$FDiv,
              'PCoA_BC' = PCoA_BC))
}
