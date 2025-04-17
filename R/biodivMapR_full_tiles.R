#' computes diversity metrics from raster
#'
#' @param dsn_grid character. path for the tiling grid
#' @param feature_dir character. path where to get features
#' @param list_features character. list of features
#' @param mask_dir character. path for masks
#' @param output_dir character. path where to save results
#' @param window_size numeric. number of clusters used in kmeans
#' @param plots list. list of sf plots
#' @param nbsamples_alpha numeric. number of samples to compute beta diversity from
#' @param nbsamples_beta numeric. number of samples to compute beta diversity from
#' @param nbCPU numeric. Number of CPUs available
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param maxRows numeric. maximum number of rows
#' @param MovingWindow boolean. should moving window be used?
#' @param siteName character. name for the output files
#'
#' @return mosaic_path
#' @export

biodivMapR_full_tiles <- function(dsn_grid, feature_dir, list_features,
                                  mask_dir = NULL, output_dir, window_size,
                                  plots, nbsamples_alpha = 1e5, nbsamples_beta = 2e3,
                                  nbCPU = 1, pcelim = 0.02, maxRows = 1200,
                                  MovingWindow = F, siteName = NULL){

  # update mask based on IQR filtering for each feature
  mask_path_list <- compute_mask_IQR_tiles(feature_dir = feature_dir,
                                           feature_list = list_features,
                                           mask_dir = mask_dir,
                                           plots = plots,
                                           nbCPU = nbCPU)

  # check which masks exist and discard plots with no masks
  mask_path <- mask_path_list$mask_path
  ID_aoi <- mask_path_list$tile_exists
  plots <- plots[ID_aoi]

  # load kmeans and beta info if exist
  Kmeans_path <- file.path(output_dir, 'Kmeans_info.RData')
  Beta_path <- file.path(output_dir, 'Beta_info.RData')
  if (file.exists(Kmeans_path)) load(Kmeans_path)
  if (file.exists(Beta_path)) load(Beta_path)
  # compute kmeans and beta info if exist
  if (!file.exists(Kmeans_path) | ! file.exists(Beta_path)){
    # define sampling points for alpha and beta diversity mapping
    samples_alpha_beta <- sample_from_plots(feature_dir = feature_dir,
                                            list_features = list_features,
                                            mask_dir = mask_dir,
                                            plots = plots,
                                            window_size = window_size,
                                            nbsamples_alpha = nbsamples_alpha,
                                            nbsamples_beta = nbsamples_beta,
                                            nbCPU = 1)
    alpha_samples <- samples_alpha_beta$samples_alpha[list_features]
    beta_samples <- samples_alpha_beta$samples_beta[c(list_features,'ID')]

    # compute alpha and beta models from samples
    nbCPU2 <- min(c(4, nbCPU))
    if (!file.exists(Kmeans_path)){
      alpha_samples$ID <- NULL
      Kmeans_info <- init_kmeans_samples(rast_sample = alpha_samples,
                                         output_dir = output_dir, nbCPU = nbCPU2)
    }
    if (!file.exists(Beta_path))
      Beta_info <- init_PCoA_samples(rast_sample = beta_samples,
                                     output_dir = output_dir,
                                     Kmeans_info = Kmeans_info,
                                     dimPCoA = 3, nbCPU = 1)
  }

  message('applying biodivMapR on tiles')
  if (nbCPU>1){
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = length(ID_aoi))
      future.apply::future_lapply(X = ID_aoi,
                                  FUN = run_biodivMapR_plot,
                                  feature_dir = feature_dir,
                                  mask_dir = mask_dir,
                                  list_features = list_features,
                                  Kmeans_info = Kmeans_info,
                                  Beta_info = Beta_info,
                                  output_dir = output_dir,
                                  window_size = window_size,
                                  maxRows = maxRows,
                                  pcelim = pcelim,
                                  MovingWindow = MovingWindow, p = p,
                                  future.seed = T,
                                  future.chunk.size = NULL,
                                  future.scheduling = structure(TRUE, ordering = "random"))
    })
    parallel::stopCluster(cl)
    plan(sequential)
  } else if (nbCPU==1){
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = length(ID_aoi))
      lapply(X = ID_aoi,
             FUN = run_biodivMapR_plot,
             feature_dir = feature_dir,
             mask_dir = mask_dir,
             list_features = list_features,
             Kmeans_info = Kmeans_info,
             Beta_info = Beta_info,
             output_dir = output_dir,
             window_size = window_size,
             maxRows = maxRows,
             MovingWindow = MovingWindow, p = p)
    })
  }

  # produce mosaic for outputs
  indices <- c('shannon', 'beta')
  mosaic_path <- list()
  for (biodividx in indices){
    # identify files
    if (biodividx == 'beta'){
      selfiles <- list.files(path = output_dir, pattern = biodividx)
    } else {
      selfiles <- list.files(path = output_dir, pattern = biodividx)
      selfiles <- selfiles[grepl(x = selfiles, pattern = "mean.tiff")]
    }
    selfiles <- file.path(output_dir, selfiles)
    # create directory
    diridx <- file.path(output_dir,biodividx)
    dir.create(diridx, showWarnings = F, recursive = T)
    # move files from - to
    files_in <- selfiles
    files_out <- as.list(file.path(diridx, basename(selfiles)))
    mapply(FUN = file.rename, from = files_in, to = files_out)
    mosaic_path[[biodividx]] <- mosaic_tiles(pattern = biodividx, siteName = siteName,
                                             dir_path = diridx, overwrite = F,
                                             vrt_save = output_dir)
  }
  return(mosaic_path)
}
