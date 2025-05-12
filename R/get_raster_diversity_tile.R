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

  res_shapeChunk <- list()
  for (idx in alphametrics){
    res_shapeChunk[[idx]] <- list()
    for (crit in c('mean', 'sd'))
      res_shapeChunk[[idx]][[crit]] <- matrix(NA,
                                              nrow = nrow(sample_newres),
                                              ncol = ncol(sample_newres))
  }
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

  # compute diversity
  if (!is.null(inputdata)){
    if (length(inputdata)>0){
      SSchunk <- get_spectralSpecies(inputdata = inputdata,
                                     Kmeans_info = Kmeans_info,
                                     selected_bands = selected_bands)
      # 5- split data chunk by window and by nbCPU to ensure parallel computing
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
      # 7- reshape alpha diversity metrics
      IDwindow <- unlist(SSwindows_per_CPU$IDwindow_perCPU)
      for (idx in alphametrics){
        for (crit in c('mean', 'sd')){
          elem <- paste0(idx,'_',crit)
          restmp <- unlist(lapply(alphabetaIdx,'[[',elem))
          res_shapeChunk[[idx]][[crit]][IDwindow] <- restmp
        }
      }
      if (!is.null(Beta_info) & !is.null(IDwindow) & length(IDwindow)>0) {
        PCoA_BC0 <- do.call(rbind,lapply(alphabetaIdx,'[[','PCoA_BC'))
        for (i in seq_len(dimPCO))
          PCoA_raster[[i]][IDwindow] <- PCoA_BC0[,i]
      }
    }
  }
  rm(alphabetaIdx)
  gc()
  ab_div_metrics <- list('richness' = res_shapeChunk$richness,
                         'shannon' = res_shapeChunk$shannon,
                         'simpson' = res_shapeChunk$simpson,
                         'fisher' = res_shapeChunk$fisher,
                         'hill' = res_shapeChunk$hill,
                         'FRic' = res_shapeChunk$FRic,
                         'FEve' = res_shapeChunk$FEve,
                         'FDiv' = res_shapeChunk$FDiv,
                         'FDis' = res_shapeChunk$FDis,
                         'FRaoq' = res_shapeChunk$FRaoq,
                         'PCoA_BC' = PCoA_raster)
  return(ab_div_metrics)
}
