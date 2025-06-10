#' Computes diversity metrics from raster data
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
#' @param maxRows numeric. max number of rows processed once by each CPU
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#'
#' @return ab_div_metrics list. contains all metrics
#' @import cli
#' @importFrom terra rast blocks readStart readStop
#' @importFrom progressr progressor handlers with_progress
#' @export

get_raster_diversity <- function(input_raster_path, Kmeans_info, Beta_info,
                                 input_mask_path = NULL, selected_bands = NULL,
                                 alphametrics = 'shannon', Hill_order = 1,
                                 FDmetric = NULL, window_size, maxRows = NULL,
                                 pcelim = 0.02, nbCPU = 1, min_sun = 0.25){
  if (is.null(maxRows))
    maxRows <- 20*window_size
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
  # Adjust size of individual chunks
  brast <- terra::rast(input_raster_path[[1]])
  blk <- terra::blocks(brast[[1]], n = 10)
  blk <- maxRows_chunk(blk = blk, maxRows = maxRows)
  blk <- nbRows_chunk(blk = blk, nbRows = window_size)
  for (fid in names(r_in))
    terra::readStart(r_in[[fid]])
  # produce a list of blocks to read and process
  blk_list <- list()
  nbchunks <- length(blk$row)
  for (i in seq_len(nbchunks))
    blk_list[[i]] <- list('row'= blk$row[i], 'nrows' = blk$nrows[i])

  if (nbCPU==1){
    # compute diversity metrics for each block
    ab_div_metrics <- lapply(X = blk_list,
                             FUN = biodivMapR_chunk,
                             Kmeans_info = Kmeans_info,
                             Beta_info = Beta_info,
                             alphametrics = alphametrics,
                             Hill_order = Hill_order,
                             FDmetric = FDmetric,
                             r_in = r_in,
                             window_size = window_size,
                             selected_bands = selected_bands,
                             pcelim = pcelim, nbCPU = nbCPU,
                             min_sun = min_sun)
  } else {
    if (nbCPU>1){
      # compute diversity metrics for each block
      # progressr::handlers(global = TRUE)
      suppressWarnings(progressr::handlers("cli"))
      # progressr::handlers("debug")
      suppressWarnings(with_progress({
        p <- progressr::progressor(steps = length(blk_list))
        ab_div_metrics <- lapply(X = blk_list,
                                 FUN = biodivMapR_chunk,
                                 Kmeans_info = Kmeans_info,
                                 Beta_info = Beta_info,
                                 alphametrics = alphametrics,
                                 Hill_order = Hill_order,
                                 FDmetric = FDmetric,
                                 r_in = r_in,
                                 window_size = window_size,
                                 selected_bands = selected_bands,
                                 pcelim = pcelim, nbCPU = nbCPU,
                                 min_sun = min_sun, p = p)
      }))
    }
  }
  for (fid in names(r_in))
    terra::readStop(r_in[[fid]])
  return(ab_div_metrics)
}
