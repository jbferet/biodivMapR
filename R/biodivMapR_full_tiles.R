#' computes diversity metrics from raster
#'
#' @param feature_dir character. path where to get features
#' @param list_features character. list of features
#' @param mask_dir character. path for masks
#' @param output_dir character. path where to save results
#' @param window_size numeric. number of clusters used in kmeans
#' @param plots list. list of sf plots
#' @param nbCPU numeric. Number of CPUs available
#' @param siteName character. name for the output files
#' @param options list. including
#' - nb_clusters numeric. number of clusters
#' - nb_samples_alpha numeric. number of samples to compute alpha diversity
#' - nb_samples_beta numeric. number of samples to compute beta diversity
#' - alphametrics character.
#' - Hill_order numeric.
#' - FDmetric character.
#' - nb_iter numeric. Number of iterations required to compute diversity
#' - pcelim numeric. minimum proportion of pixels to consider spectral species
#' - maxRows numeric. maximum number of rows
#' - moving_window boolean. should moving window be used?
#' - mosaic_output boolean. set TRUE if outputs need to be mosaiced
#' - weightIRQ numeric. IQR applied to filter out features to be used
#'
#' @return mosaic_path
#' @export

biodivMapR_full_tiles <- function(feature_dir, list_features, mask_dir = NULL,
                                  output_dir, window_size, plots, nbCPU = 1,
                                  siteName = NULL, options = NULL){

  # define options
  options <- set_options_biodivMapR(fun = 'biodivMapR_full_tiles', options = options)
  nb_clusters <- options$nb_clusters
  nb_samples_alpha<- options$nb_samples_alpha
  nb_samples_beta <- options$nb_samples_beta
  alphametrics <- options$alphametrics
  Hill_order <- options$Hill_order
  FDmetric <- options$FDmetric
  nb_iter <- options$nb_iter
  pcelim <- options$pcelim
  maxRows <- options$maxRows
  moving_window <- options$moving_window
  mosaic_output <- options$mosaic_output
  weightIRQ <- options$weightIRQ

  # # adjust hill name if required
  # if ('hill' %in% alphametrics){
  #
  # }

  # sample data if not already sampled
  samples <- biodivMapR_sample(feature_dir = feature_dir,
                               list_features = list_features,
                               mask_dir = mask_dir,
                               output_dir = output_dir,
                               window_size = window_size,
                               plots = plots,
                               nb_clusters = nb_clusters,
                               nb_samples_alpha = nb_samples_alpha,
                               nb_samples_beta = nb_samples_beta,
                               pcelim = pcelim, nbCPU = nbCPU,
                               nb_iter = nb_iter,
                               weightIRQ = weightIRQ)

  message('applying biodivMapR on tiles')
  maxCPU <- length(samples$ID_aoi)
  if (nbCPU > maxCPU)
    nbCPU <-  maxCPU

  if (nbCPU>1){
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = maxCPU)
      future.apply::future_lapply(X = samples$ID_aoi,
                                  FUN = run_biodivMapR_plot,
                                  feature_dir = feature_dir,
                                  mask_dir = mask_dir,
                                  list_features = list_features,
                                  Kmeans_info = samples$Kmeans_info,
                                  Beta_info = samples$Beta_info,
                                  alphametrics = alphametrics,
                                  Hill_order = Hill_order,
                                  FDmetric = FDmetric,
                                  output_dir = output_dir,
                                  window_size = window_size,
                                  maxRows = maxRows,
                                  pcelim = pcelim,
                                  moving_window = moving_window, p = p,
                                  future.seed = TRUE,
                                  future.chunk.size = NULL,
                                  future.scheduling = structure(TRUE,
                                                                ordering = "random"))
    })
    parallel::stopCluster(cl)
    plan(sequential)
  } else if (nbCPU==1){
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = maxCPU)
      lapply(X = samples$ID_aoi,
             FUN = run_biodivMapR_plot,
             feature_dir = feature_dir,
             mask_dir = mask_dir,
             list_features = list_features,
             Kmeans_info = samples$Kmeans_info,
             Beta_info = samples$Beta_info,
             alphametrics = alphametrics,
             Hill_order = Hill_order,
             FDmetric = FDmetric,
             output_dir = output_dir,
             window_size = window_size,
             maxRows = maxRows,
             moving_window = moving_window, p = p)
    })
  }

  if (mosaic_output){
    # produce mosaic for outputs
    indices <- c(alphametrics, 'beta', FDmetric)
    mosaic_path <- list()
    for (biodividx in indices){
      # create directory
      diridx <- file.path(output_dir,biodividx)
      dir.create(diridx, showWarnings = FALSE, recursive = TRUE)
      # identify files
      selfiles <- list.files(path = output_dir, pattern = biodividx)
      if (! biodividx %in% alphametrics){
        selfiles <- file.path(output_dir, selfiles)
        # move files from - to
        files_in <- selfiles
        files_out <- as.list(file.path(diridx, basename(selfiles)))
        mapply(FUN = file.rename, from = files_in, to = files_out)
        mosaic_path[[biodividx]] <- mosaic_tiles(pattern = biodividx,
                                                 siteName = siteName,
                                                 dir_path = diridx,
                                                 overwrite = FALSE,
                                                 vrt_save = output_dir)
      } else if (biodividx %in% alphametrics){
        # compute mean mosaic
        selfiles_mean <- selfiles[grepl(x = basename(selfiles),
                                        pattern = "mean.tiff")]
        selfiles_mean <- file.path(output_dir, selfiles_mean)
        # move files from - to
        files_in <- selfiles_mean
        files_out <- as.list(file.path(diridx, basename(selfiles_mean)))
        mapply(FUN = file.rename, from = files_in, to = files_out)
        mosaic_path[[biodividx]] <- mosaic_tiles(pattern = 'mean.tiff',
                                                 siteName = siteName,
                                                 dir_path = diridx,
                                                 overwrite = FALSE,
                                                 vrt_save = output_dir)

        # compute sd mosaic
        selfiles_sd <- selfiles[grepl(x = basename(selfiles),
                                        pattern = 'sd.tiff')]
        selfiles_sd <- file.path(output_dir, selfiles_sd)
        # move files from - to
        files_in <- selfiles_sd
        files_out <- as.list(file.path(diridx, basename(selfiles_sd)))
        mapply(FUN = file.rename, from = files_in, to = files_out)
        mosaic_path[[biodividx]] <- mosaic_tiles(pattern = 'sd.tiff',
                                                 siteName = siteName,
                                                 dir_path = diridx,
                                                 overwrite = FALSE,
                                                 vrt_save = output_dir)
      }
    }
  }
  return(mosaic_path)
}
