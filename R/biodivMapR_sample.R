#' sample pixels and plots to produce alpha and beta models, including
#' clustering, dissimilarity & PCoA
#'
#' @param feature_dir character. path where to get features
#' @param list_features character. list of features
#' @param mask_dir character. path for masks
#' @param output_dir character. path where to save results
#' @param window_size numeric. number of clusters used in kmeans
#' @param plots list. list of sf plots
#' @param nb_clusters numeric. number of clusters
#' @param nb_samples_alpha numeric. number of samples to compute alpha diversity
#' @param beta_metrics boolean. set TRUE to compute beta diversity
#' @param nb_samples_beta numeric. number of samples to compute beta diversity
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param nb_iter numeric. Number of iterations required to compute diversity
#' @param weightIRQ numeric. IQR applied to filter out features to be used
#'
#' @return mosaic_path
#' @export

biodivMapR_sample <- function(feature_dir, list_features, mask_dir = NULL,
                              output_dir, window_size, plots, nb_clusters = 50,
                              nb_samples_alpha = 1e5, beta_metrics = TRUE,
                              nb_samples_beta = 2e3, pcelim = 0.02, nbCPU = 1,
                              nb_iter = 10, weightIRQ = 4){

  message('biodivMapR sampling')
  # update mask based on IQR filtering for each feature
  mask_path_list <- compute_mask_iqr_tiles(feature_dir = feature_dir,
                                           feature_list = list_features,
                                           mask_dir = mask_dir,
                                           plots = plots,
                                           nbCPU = nbCPU,
                                           weightIRQ = weightIRQ)
  gc()
  # check which masks exist and discard plots with no masks
  ID_aoi <- mask_path_list$tile_exists
  plots <- plots[ID_aoi]

  # load kmeans and beta info if exist
  Kmeans_path <- file.path(output_dir, 'Kmeans_info.RData')
  Beta_path <- file.path(output_dir, 'Beta_info.RData')
  Beta_info <- NULL
  if (file.exists(Kmeans_path))
    load(Kmeans_path)
  if (file.exists(Beta_path) & beta_metrics)
    load(Beta_path)
  # compute kmeans and beta info if exist
  if (!file.exists(Kmeans_path) | (! file.exists(Beta_path) & beta_metrics)){
    # define sampling points for alpha and beta diversity mapping
    samples_alpha_beta <- sample_from_plots(feature_dir = feature_dir,
                                            list_features = list_features,
                                            mask_dir = mask_dir,
                                            plots = plots,
                                            window_size = window_size,
                                            nb_samples_alpha = nb_samples_alpha,
                                            beta_metrics = beta_metrics,
                                            nb_samples_beta = nb_samples_beta,
                                            nbCPU = nbCPU)
    alpha_samples <- samples_alpha_beta$samples_alpha[list_features]
    if (beta_metrics)
      beta_samples <- samples_alpha_beta$samples_beta[c(list_features,'ID')]

    # compute alpha and beta models from samples
    nbCPU2 <- min(c(4, nbCPU))
    if (!file.exists(Kmeans_path)){
      alpha_samples$ID <- NULL
      Kmeans_info <- init_kmeans_samples(rast_sample = alpha_samples,
                                         output_dir = output_dir,
                                         nb_clusters = nb_clusters,
                                         nb_iter = nb_iter,
                                         nbCPU = nbCPU2)
    }
    if (!file.exists(Beta_path) & beta_metrics)
      Beta_info <- init_PCoA_samples(rast_sample = beta_samples,
                                     output_dir = output_dir,
                                     Kmeans_info = Kmeans_info,
                                     pcelim = pcelim,
                                     dimPCoA = 3, nbCPU = 1)
  }
  message('sampling succeeded')
  return(list('Kmeans_info' = Kmeans_info, 'Beta_info' = Beta_info,
              'ID_aoi' = ID_aoi))
}
