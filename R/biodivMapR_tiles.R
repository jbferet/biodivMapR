#' computes diversity metrics from raster
#'
#' @param feature_dir character. path where to get features
#' @param list_features character. list of features
#' @param mask_dir character. path for masks
#' @param output_dir character. path where to save results
#' @param window_size numeric. number of clusters used in kmeans
#' @param plots list. list of sf plots
#' @param alpha_metrics character.
#' @param Hill_order numeric.
#' @param fd_metrics character.
#' @param nbCPU numeric. Number of CPUs available
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param maxRows numeric. maximum number of rows
#' @param moving_window boolean. should moving window be used?
#' @param site_name character. name for the output files
#' @param mosaic_output boolean. set TRUE if outputs need to be mosaiced
#'
#' @return mosaic_path
#' @importFrom parallel makeCluster stopCluster
#' @importFrom future plan sequential
#' @importFrom future.apply future_mapply
#' @importFrom progressr with_progress progressor handlers
#' @export

biodivMapR_tiles <- function(feature_dir, list_features, mask_dir = NULL,
                             output_dir, window_size, plots,
                             alpha_metrics = 'shannon', Hill_order = 1,
                             fd_metrics = NULL, nbCPU = 1, pcelim = 0.02,
                             maxRows = 1200, moving_window = FALSE,
                             site_name = NULL, mosaic_output = TRUE){

  # update mask based on IQR filtering for each feature
  mask_path_list <- compute_mask_iqr_tiles(feature_dir = feature_dir,
                                           feature_list = list_features,
                                           mask_dir = mask_dir,
                                           plots = plots,
                                           nbCPU = nbCPU)

  # check which masks exist and discard plots with no masks
  ID_aoi <- mask_path_list$tile_exists
  plots <- plots[ID_aoi]

  # load kmeans and beta info if exist
  Kmeans_path <- file.path(output_dir, 'Kmeans_info.RData')
  Beta_path <- file.path(output_dir, 'Beta_info.RData')
  Kmeans_info <- Beta_info <- NULL
  if (file.exists(Kmeans_path) & file.exists(Beta_path)){
    load(Kmeans_path)
    load(Beta_path)
  } else {
    message('please perform sampling with function "biodivMapR_sample"')
    message('"Kmeans_info.RData" & "Beta_info.RData" missing. stopping process')
    stop()
  }

  message('applying biodivMapR on tiles')
  maxCPU <- length(ID_aoi)
  if (nbCPU > maxCPU)
    nbCPU <-  maxCPU

  if (nbCPU>1){
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = maxCPU)
      future.apply::future_lapply(X = ID_aoi, FUN = run_biodivMapR_plot,
                                  feature_dir = feature_dir,
                                  mask_dir = mask_dir,
                                  list_features = list_features,
                                  Kmeans_info = Kmeans_info,
                                  Beta_info = Beta_info,
                                  alpha_metrics = alpha_metrics,
                                  Hill_order = Hill_order,
                                  fd_metrics = fd_metrics, output_dir = output_dir,
                                  window_size = window_size,
                                  maxRows = maxRows, pcelim = pcelim,
                                  moving_window = moving_window, p = p,
                                  future.seed = TRUE, future.chunk.size = NULL,
                                  future.scheduling = structure(TRUE,
                                                                ordering = "random"))
    })
    parallel::stopCluster(cl)
    plan(sequential)
  } else if (nbCPU==1){
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = maxCPU)
      lapply(X = ID_aoi, FUN = run_biodivMapR_plot,
             feature_dir = feature_dir, mask_dir = mask_dir,
             list_features = list_features, Kmeans_info = Kmeans_info,
             Beta_info = Beta_info, alpha_metrics = alpha_metrics,
             Hill_order = Hill_order, fd_metrics = fd_metrics,
             output_dir = output_dir, window_size = window_size,
             maxRows = maxRows, moving_window = moving_window, p = p)
    })
  }

  if (mosaic_output){
    # produce mosaic for outputs
    indices <- c(alpha_metrics, 'beta', fd_metrics)
    mosaic_path <- list()
    for (biodividx in indices){
      # identify files
      if (! biodividx %in% alpha_metrics){
        selfiles <- list.files(path = output_dir, pattern = biodividx)
      } else if (biodividx %in% alpha_metrics){
        selfiles <- list.files(path = output_dir, pattern = biodividx)
        selfiles <- selfiles[grepl(x = basename(selfiles),
                                   pattern = "mean.tiff")]
      }
      selfiles <- file.path(output_dir, selfiles)
      # create directory
      diridx <- file.path(output_dir,biodividx)
      dir.create(diridx, showWarnings = FALSE, recursive = TRUE)
      # move files from - to
      files_in <- selfiles
      files_out <- as.list(file.path(diridx, basename(selfiles)))
      mapply(FUN = file.rename, from = files_in, to = files_out)
      mosaic_path[[biodividx]] <- mosaic_tiles(pattern = biodividx,
                                               site_name = site_name,
                                               dir_path = diridx,
                                               overwrite = FALSE,
                                               vrt_save = output_dir)
    }
  }
  return(mosaic_path)
}
