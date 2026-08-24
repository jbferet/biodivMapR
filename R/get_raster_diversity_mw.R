#' Computes diversity metrics from raster data based on moving window
#'
#' @param input_raster_path list. list of paths corresponding to input rasters
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param input_mask_path character. path for mask file
#' @param selected_bands numeric. bands selected from input_rast
#' @param alpha_metrics list. alpha diversity metrics: richness, shannon, simpson
#' @param Hill_order numeric. Hill order
#' @param beta_metrics list. beta diversity metrics: bray, brayturn,
#' @param fd_metrics character. list of functional metrics
#' @param window_size numeric. window size for square plots
#' @param maxRows numeric. max number of rows in each block
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#'
#' @return ab_div_metrics list. contains all metrics
#' @importFrom terra rast
#' @export

get_raster_diversity_mw <- function(
    input_raster_path, Kmeans_info, Beta_info, input_mask_path = NULL,
    selected_bands = NULL, alpha_metrics = 'shannon', Hill_order = 1,
    beta_metrics = 'bray', fd_metrics = NULL, window_size, maxRows = NULL,
    pcelim = 0.02, nbCPU = 1, min_sun = 0.25){

  if (!is.null(fd_metrics)){
    message('functional diversity metrics are not too computationally intensive for moving windows')
    message('setting fd_metrics = FALSE')
    fd_metrics <- FALSE
  }

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

  res_shapeChunk <- list()
  for (idx in alpha_metrics){
    res_shapeChunk[[idx]] <- list()
    for (crit in c('mean', 'sd'))
      res_shapeChunk[[idx]][[crit]] <- matrix(NA,
                                              nrow = nrow(rast_in),
                                              ncol = ncol(rast_in))
  }
  dimPCO <- 3
  if (!is.null(Beta_info))
    dimPCO <- ncol(Beta_info$beta_pco[[1]]$points)
  pcoa_raster <- list()
  for (beta in beta_metrics)
    pcoa_raster[[beta]] <- list('PCoA1' = matrix(NA, nrow = nrow(rast_in),
                                                 ncol = ncol(rast_in)),
                                'PCoA2' = matrix(NA, nrow = nrow(rast_in),
                                                 ncol = ncol(rast_in)),
                                'PCoA3' = matrix(NA, nrow = nrow(rast_in),
                                                 ncol = ncol(rast_in)))

  # v1: line per line
  nbWindows <- ncol(rast_in)
  rastmat <- as.matrix(rast_in[[1]], wide = TRUE)
  whichj <- which(!rowSums(rastmat, na.rm = TRUE)==0)
  for (j in whichj){
    inputdata <- list()
    sel <- NULL
    notNA <- which(!is.na(rast_in[j,,1]))
    if (length(notNA)>0){
      for (i in notNA){
        xmin <- max(c(1, i-(window_size-1)/2))
        xmax <- min(c(ncol(rast_in), i+(window_size-1)/2))
        ymin <- max(c(1, j-(window_size-1)/2))
        ymax <- min(c(nrow(rast_in), j+(window_size-1)/2))
        inputdata[[i]] <- na.omit(rast_in[ymin:ymax, xmin:xmax, ])
        if (!is.null(inputdata[[i]]$mask)){
          sel <- which(inputdata[[i]]$mask==1)
          if (length(sel)>0)
            inputdata[[i]] <- inputdata[[i]][sel,]
          inputdata[[i]]$mask <- NULL
        }
        if (nrow(inputdata[[i]])>min_sun*window_size**2){
          inputdata[[i]]$win_ID <- i
        } else {
          inputdata[[i]] <- NULL
        }
      }
      inputdata <- do.call(what = rbind, args = inputdata)
    }

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
        alphabetaIdx <- lapply(X = SSwindows_per_CPU$SSwindow_perCPU,
                               # FUN = alphabeta_window_list,
                               FUN = alphabeta_chunk,
                               nb_clusters = nb_clusters,
                               alpha_metrics = alpha_metrics,
                               Hill_order = Hill_order,
                               Beta_info = Beta_info,
                               beta_metrics = beta_metrics,
                               pcelim = pcelim)

        # 7- reshape alpha diversity metrics
        IDwindow <- unlist(SSwindows_per_CPU$IDwindow_perCPU)
        for (idx in alpha_metrics){
          for (crit in c('mean', 'sd')){
            elem <- paste0(idx,'_',crit)
            restmp <- unlist(lapply(alphabetaIdx,'[[',elem))
            res_shapeChunk_tmp <- rep(x = NA, nbWindows)
            if (length(IDwindow)>0)
              res_shapeChunk_tmp[IDwindow] <- restmp
            res_shapeChunk[[idx]][[crit]][j,] <- res_shapeChunk_tmp
          }
        }

        # resahape beta
        for (beta in beta_metrics){
          pcoa_diss <- matrix(data = NA, nrow = nbWindows, ncol = dimPCO)
          if (!is.null(Beta_info) & !is.null(IDwindow) & length(IDwindow)>0) {
            pcoa_diss0 <- do.call(rbind,lapply(alphabetaIdx,'[[',
                                               paste0('pcoa_', beta)))
            pcoa_diss[IDwindow,] <- pcoa_diss0
            for (i in 1:dimPCO)
              pcoa_raster[[beta]][[i]][j,] <- pcoa_diss[,i]
          }
        }
      }
    }
  }


  rm(alphabetaIdx)
  gc()
  ab_div_metrics <- list('richness' = res_shapeChunk$richness,
                         'shannon' = res_shapeChunk$shannon,
                         'simpson' = res_shapeChunk$simpson,
                         'hill' = res_shapeChunk$hill,
                         'FRic' = res_shapeChunk$FRic,
                         'FEve' = res_shapeChunk$FEve,
                         'FDiv' = res_shapeChunk$FDiv,
                         'FDis' = res_shapeChunk$FDis,
                         'FRaoq' = res_shapeChunk$FRaoq,
                         'pcoa_bray' = pcoa_raster[['bray']],
                         'pcoa_brayturn' = pcoa_raster[['brayturn']],
                         'pcoa_simpson_diss' = pcoa_raster[['simpson_diss']],
                         'pcoa_jaccard' = pcoa_raster[['jaccard']],
                         'pcoa_jaccardturn' = pcoa_raster[['jaccardturn']],
                         'pcoa_sorensen' = pcoa_raster[['sorensen']])
  return(ab_div_metrics)
}
