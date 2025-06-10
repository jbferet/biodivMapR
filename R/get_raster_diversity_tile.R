#' Computes diversity metrics from raster data based on moving window
#'
#' @param input_raster_path list. list of paths corresponding to input rasters
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param input_mask_path character. path for mask file
#' @param selected_bands numeric. bands selected from input_rast
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param Hill_order numeric. Hill order
#' @param FDmetric character. list of functional metrics
#' @param window_size numeric. window size for square plots
#' @param maxRows numeric. max number of rows in each block
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#'
#' @return ab_div_metrics list. contains all metrics
#' @import cli
#' @importFrom terra rast blocks readStart readStop
#' @importFrom progressr progressor handlers with_progress
#' @export

get_raster_diversity_tile <- function(input_raster_path, Kmeans_info, Beta_info,
                                      input_mask_path = NULL, selected_bands = NULL,
                                      alphametrics = 'shannon', Hill_order = 1,
                                      FDmetric = NULL, window_size, maxRows = NULL,
                                      pcelim = 0.02, nbCPU = 1, min_sun = 0.25){
  # prepare to read input raster data
  r_in <- list()
  if (is.null(names(input_raster_path)))
    names(input_raster_path) <- seq_len(length(input_raster_path))
  for (fid in names(input_raster_path))
    r_in[[fid]] <- terra::rast(input_raster_path[[fid]])
  # if a mask file is provided
  if (!is.null(input_mask_path)) {
    r_in[['mask']] <- terra::rast(input_mask_path)
    names(r_in[['mask']]) <- 'mask'
    input_raster_path[['mask']] <- input_mask_path
  }
  rast_in <- lapply(X = input_raster_path, FUN = terra::rast)
  rast_in <- terra::rast(rast_in)
  if (is.null(input_mask_path)) {
    selected_bands <- seq_len(dim(rast_in)[3])
  } else {
    selected_bands <- seq_len(dim(rast_in)[3]-1)
  }

  # define template for diversity rasters
  sample_newres <- terra::aggregate(x = rast_in[[1]], fact = window_size)
  terra::values(sample_newres) <- seq_len(length(terra::values(sample_newres)))
  sample_tmp <- terra::resample(x = sample_newres,
                                y = rast_in[[1]],
                                method = 'near')
  names(sample_tmp) <- 'win_ID'

  res_shape_chunk <- list()
  for (idx in alphametrics){
    res_shape_chunk[[idx]] <- list()
    for (crit in c('mean', 'sd'))
      res_shape_chunk[[idx]][[crit]] <- matrix(NA,
                                               nrow = nrow(sample_newres),
                                               ncol = ncol(sample_newres))
  }
  for (idx in FDmetric)
    res_shape_chunk[[idx]] <- matrix(NA,
                                     nrow = nrow(sample_newres),
                                     ncol = ncol(sample_newres))
  dimPCO <- 3
  if (!is.null(Beta_info))
    dimPCO <- ncol(Beta_info$BetaPCO$points)
  PCoA_raster <- list('PCoA1' = matrix(NA, nrow = nrow(sample_newres),
                                       ncol = ncol(sample_newres)),
                      'PCoA2' = matrix(NA, nrow = nrow(sample_newres),
                                       ncol = ncol(sample_newres)),
                      'PCoA3' = matrix(NA, nrow = nrow(sample_newres),
                                       ncol = ncol(sample_newres)))


  rast_in <- c(rast_in, sample_tmp)
  # v1: line per line
  inputdata <- as.data.frame(terra::values(rast_in))
  names(inputdata) <- names(rast_in)
  inputdata <- na.omit(inputdata)
  abund <- table(inputdata$win_ID)
  selWin <- as.numeric(names(abund)[which(abund > min_sun*window_size**2)])
  selpix <- which(inputdata$win_ID %in% selWin)
  inputdata <- inputdata[selpix, ]
  nbWindows <- max(inputdata$win_ID)

  # compute diversity
  if (!is.null(inputdata)){
    if (length(inputdata)>0){
      SSchunk <- get_spectralSpecies(inputdata = inputdata,
                                     Kmeans_info = Kmeans_info,
                                     selected_bands = selected_bands)
      # 5- split data chunk by window and by nbCPU to ensure parallel computing
      if (dim(SSchunk)[1]>0){
        SSwindows_per_CPU <- split_chunk(SSchunk, nbCPU)
        # 6- compute diversity metrics
        nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
        alphabetaIdx_CPU <- lapply(X = SSwindows_per_CPU$SSwindow_perCPU,
                                   FUN = alphabeta_window_list,
                                   nb_clusters = nb_clusters,
                                   alphametrics = alphametrics,
                                   Beta_info = Beta_info,
                                   Hill_order = Hill_order,
                                   pcelim = pcelim)
        alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = FALSE)


        if (!is.null(FDmetric)){
          inputdata <- cbind(center_reduce(x = inputdata[selected_bands],
                                           m = Kmeans_info$MinVal,
                                           sig = Kmeans_info$Range),
                             'win_ID' = inputdata$win_ID)
          windows_per_CPU <- split_chunk(inputdata, nbCPU)
          funct_idx_cpu <- future.apply::future_lapply(X = windows_per_CPU$SSwindow_perCPU,
                                                       FUN = functional_window_list,
                                                       FDmetric = FDmetric,
                                                       future.seed = TRUE)
          FunctionalIdx <- unlist(funct_idx_cpu, recursive = FALSE)
        }


        # 7- reshape alpha diversity metrics
        IDwindow <- unlist(SSwindows_per_CPU$IDwindow_perCPU)
        for (idx in alphametrics){
          for (crit in c('mean', 'sd')){
            elem <- paste0(idx,'_',crit)
            restmp <- unlist(lapply(alphabetaIdx,'[[',elem))
            res_shape_chunk[[idx]][[crit]][IDwindow] <- restmp
          }
        }
        if (!is.null(Beta_info) & !is.null(IDwindow) & length(IDwindow)>0) {
          PCoA_BC0 <- do.call(rbind,lapply(alphabetaIdx,'[[','PCoA_BC'))
          for (i in seq_len(dimPCO))
            PCoA_raster[[i]][IDwindow] <- PCoA_BC0[,i]
        }

        list_funct_idx <- c('FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')
        for (idx in list_funct_idx){
          if (idx %in% FDmetric){
            restmp <- unlist(lapply(FunctionalIdx,'[[',idx))
            res_shape_chunk[[idx]][IDwindow] <- restmp
            # res_shape_chunk_tmp <- rep(x = NA,nbWindows)
            # res_shape_chunk_tmp[IDwindow] <- restmp
            # res_shape_chunk[[idx]] <- matrix(res_shape_chunk_tmp,
            #                                  nrow = ceiling(blk$nrows/window_size),
            #                                  byrow = TRUE)
          }
        }
      }
    }
  }
  rm(alphabetaIdx)
  gc()
  ab_div_metrics <- list('richness' = res_shape_chunk$richness,
                         'shannon' = res_shape_chunk$shannon,
                         'simpson' = res_shape_chunk$simpson,
                         'fisher' = res_shape_chunk$fisher,
                         'hill' = res_shape_chunk$hill,
                         'FRic' = res_shape_chunk$FRic,
                         'FEve' = res_shape_chunk$FEve,
                         'FDiv' = res_shape_chunk$FDiv,
                         'FDis' = res_shape_chunk$FDis,
                         'FRaoq' = res_shape_chunk$FRaoq,
                         'PCoA_BC' = PCoA_raster)
  return(ab_div_metrics)
}
