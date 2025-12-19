#' set options
#'
#' @param fun character. name of the function which has optional parameters
#' @param options list. including
#' - nb_clusters numeric. number of clusters
#' - nb_samples_alpha numeric. number of samples to compute alpha diversity
#' - nb_samples_beta numeric. number of samples to compute beta diversity
#' - alpha_metrics character.
#' - Hill_order numeric.
#' - beta_metrics boolean.
#' - fd_metrics character.
#' - nb_iter numeric. Number of iterations required to compute diversity
#' - pcelim numeric. minimum proportion of pixels to consider spectral species
#' - maxRows numeric. maximum number of rows
#' - moving_window boolean. should moving window be used?
#' - mosaic_output boolean. set TRUE if outputs need to be mosaiced
#' - weightIQR numeric. IQR applied to filter out features to be used
#'
#' @return options with default values when missing
#' @export

set_options_biodivMapR <- function(fun, options = NULL){

  if (fun == 'biodivMapR_full'){
    if (is.null(options$alpha_metrics))
      options$alpha_metrics <- 'shannon'
    if (is.null(options$Hill_order))
      options$Hill_order <- 1
    if (is.null(options$beta_metrics))
      options$beta_metrics <- TRUE
    if (is.null(options$fd_metrics))
      options$fd_metrics <- NULL
    if (is.null(options$nb_samples_alpha))
      options$nb_samples_alpha <- 1e5
    if (is.null(options$nb_samples_beta))
      options$nb_samples_beta <- 2e3
    if (is.null(options$nb_clusters))
      options$nb_clusters <- 50
    if (is.null(options$nb_iter))
      options$nb_iter <- 10
    if (is.null(options$pcelim))
      options$pcelim <- 0.02
    if (is.null(options$maxRows))
      options$maxRows <- 1200
    if (is.null(options$moving_window))
      options$moving_window <- FALSE
    if (is.null(options$min_sun))
      options$min_sun <- 0.25
    if (is.null(options$dimPCoA))
      options$dimPCoA <- 3
    if (is.null(options$progressbar))
      options$progressbar <- TRUE
    if (is.null(options$filetype))
      options$filetype <- 'GTiff'
  }

  if (fun == 'biodivMapR_full_tiles'){
    if (is.null(options$nb_clusters))
      options$nb_clusters <- 50
    if (is.null(options$nb_samples_alpha))
      options$nb_samples_alpha <- 1e5
    if (is.null(options$nb_samples_beta))
      options$nb_samples_beta <- 2e3
    if (is.null(options$alpha_metrics))
      options$alpha_metrics <- 'shannon'
    if (is.null(options$Hill_order))
      options$Hill_order <- 1
    if (is.null(options$beta_metrics))
      options$beta_metrics <- TRUE
    if (is.null(options$fd_metrics))
      options$fd_metrics <- NULL
    if (is.null(options$nb_iter))
      options$nb_iter <- 10
    if (is.null(options$pcelim))
      options$pcelim <- 0.02
    if (is.null(options$maxRows))
      options$maxRows <- 1200
    if (is.null(options$moving_window))
      options$moving_window <- FALSE
    if (is.null(options$mosaic_output))
      options$mosaic_output <- TRUE
    if (is.null(options$weightIQR))
      options$weightIQR <- 4
  }
  return(options)
}
