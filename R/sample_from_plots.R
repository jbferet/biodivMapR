#' get samples for alpha and beta diversity mapping
#'
#' @param feature_dir character.
#' @param list_features character.
#' @param plots list.
#' @param mask_dir character.
#' @param window_size numeric.
#' @param nbCPU numeric.
#' @param nb_samples_alpha numeric.
#' @param nb_samples_beta numeric.
#'
#' @return samples_alpha_beta
#' @export

sample_from_plots <- function(feature_dir, list_features, plots, mask_dir = NULL,
                              window_size, nbCPU = 1, nb_samples_alpha = 1e5,
                              nb_samples_beta = 2e3){

  nb_pix_valid <- get_valid_pixels_from_tiles(feature_dir, plots, nbCPU = 1,
                                            mask_dir = mask_dir)
  # sample the plot network for alpha diversity
  samples_alpha_terra <- sample_from_plots_alpha(feature_dir = feature_dir,
                                                 list_features = list_features,
                                                 plots = plots,
                                                 nb_pix_valid = nb_pix_valid,
                                                 mask_dir = mask_dir,
                                                 nb_samples_alpha = nb_samples_alpha,
                                                 nbCPU = nbCPU)

  # sample the plot network for beta diversity
  samples_beta_terra <- sample_from_plots_beta(feature_dir = feature_dir,
                                               list_features = list_features,
                                               plots = plots,
                                               window_size = window_size,
                                               nb_pix_valid = nb_pix_valid,
                                               mask_dir = mask_dir,
                                               nb_samples_beta = nb_samples_beta,
                                               nbCPU = nbCPU)

  # use a unique ID per plot
  idx <- 0
  for (i in names(plots)){
    if (!is.null(samples_beta_terra[[i]])){
      for (j in seq_len(length(samples_beta_terra[[i]]))){
        if (!is.null(samples_beta_terra[[i]][[j]])){
          idx <- idx + 1
          samples_beta_terra[[i]][[j]]$ID <- idx
        }
      }
      samples_beta_terra[[i]] <- do.call('rbind', samples_beta_terra[[i]])
    }
  }
  samples_beta_terra <- do.call('rbind', samples_beta_terra)
  samples_alpha_beta <- list('samples_alpha' = samples_alpha_terra,
                             'samples_beta' = samples_beta_terra)
  return(samples_alpha_beta)
}
