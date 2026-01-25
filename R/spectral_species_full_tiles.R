#' computes diversity metrics from raster
#'
#' @param feature_dir character. path where to get features
#' @param list_features character. list of features
#' @param mask_dir character. path for masks
#' @param output_dir character. path where to save results
#' @param plots list. list of sf plots
#' @param nbCPU numeric. Number of CPUs available
#' @param site_name character. name for the output files
#' @param overwrite boolean.
#' @param options list. including
#' - nb_clusters numeric. number of clusters
#' - nb_samples_alpha numeric. number of samples to compute alpha diversity
#' - mosaic_output boolean. set TRUE if outputs need to be mosaiced
#' - weightIQR numeric. IQR applied to filter out features to be used
#'
#' @return mosaic_path
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom future plan sequential
#' @importFrom future.apply future_mapply
#' @importFrom progressr with_progress progressor handlers
#' @export

spectral_species_full_tiles <- function(feature_dir, list_features,
                                        mask_dir = NULL, output_dir,
                                        plots, nbCPU = 1, site_name = NULL,
                                        overwrite = TRUE, options = NULL){
  nb_iter <- 1
  # define options
  options <- set_options_biodivMapR(fun = 'spectral_species_full_tiles',
                                    options = options)
  nb_clusters <- options$nb_clusters
  nb_samples_alpha<- options$nb_samples_alpha
  maxRows <- options$maxRows
  weightIQR <- options$weightIQR
  Kmeans_path <- options$Kmeans_path

  # compute mask
  mask_path_list <- compute_mask_iqr_tiles(feature_dir = feature_dir,
                                           feature_list = list_features,
                                           mask_dir = mask_dir,
                                           plots = plots,
                                           nbCPU = nbCPU,
                                           weightIQR = weightIQR)

  # adjust number of clusters if less than number of plots
  maxCPU <- length(plots)
  if (nbCPU > maxCPU)
    nbCPU <-  maxCPU
  if (nbCPU > parallel::detectCores(logical = FALSE))
    nbCPU <- parallel::detectCores(logical = FALSE)

  message('alpha sampling')
  gc()
  # check which masks exist and discard plots with no masks
  ID_aoi <- mask_path_list$tile_exists
  plots <- plots[ID_aoi]
  # set default path for Kmeans_info and Beta_info
  if (is.null(Kmeans_path))
    Kmeans_path <- file.path(output_dir, 'Kmeans_info.RData')
  Kmeans_info <- NULL
  # load kmeans and beta info if exist
  if (file.exists(Kmeans_path))
    load(Kmeans_path)
  # compute kmeans and beta info if exist
  if (!file.exists(Kmeans_path)){
    # define sampling points for alpha and beta diversity mapping
    nb_pix_valid <- get_valid_pixels_from_tiles(feature_dir, plots, nbCPU = 1,
                                                mask_dir = mask_dir)
    # sample the plot network for alpha diversity
    samples_alpha <- sample_from_plots_alpha(feature_dir = feature_dir,
                                             list_features = list_features,
                                             plots = plots,
                                             nb_pix_valid = nb_pix_valid,
                                             mask_dir = mask_dir,
                                             nb_samples_alpha = nb_samples_alpha,
                                             nbCPU = nbCPU)
    # compute alpha and beta models from samples
    alpha_samples <- samples_alpha[list_features]
    alpha_samples$ID <- NULL
    Kmeans_info <- init_kmeans_samples(rast_sample = alpha_samples,
                                       output_dir = output_dir,
                                       nb_clusters = nb_clusters,
                                       nb_iter = nb_iter,
                                       nbCPU = 1)
  }
  message('applying kmeans on tiles')
  maxCPU <- length(ID_aoi)
  if (nbCPU > maxCPU)
    nbCPU <-  maxCPU

  if (nbCPU>1){
    cl <- parallel::makeCluster(nbCPU)
    with(future::plan("cluster", workers = cl), local = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = maxCPU)
      ss_path <- future.apply::future_lapply(X = ID_aoi,
                                             FUN = apply_spectral_species_plot,
                                             feature_dir = feature_dir,
                                             mask_dir = mask_dir,
                                             list_features = list_features,
                                             Kmeans_info = Kmeans_info,
                                             output_dir = output_dir,
                                             overwrite = overwrite, p = p,
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
      ss_path <- lapply(X = ID_aoi,
                        FUN = apply_spectral_species_plot, feature_dir = feature_dir,
                        mask_dir = mask_dir, list_features = list_features,
                        Kmeans_info = Kmeans_info, output_dir = output_dir,
                        overwrite = overwrite, p = p)
    })
  }
  return(ss_path)
}
