#' get samples for alpha and beta diversity mapping
#'
#' @param feature_dir character.
#' @param list_features character.
#' @param plots list.
#' @param mask_dir character.
#' @param window_size numeric.
#' @param nbCPU numeric.
#' @param nbsamples_alpha numeric.
#' @param nbsamples_beta numeric.
#'
#' @return samples_alpha_beta
#' @export

sample_from_plots <- function(feature_dir, list_features, plots, mask_dir = NULL,
                              window_size, nbCPU = 1, nbsamples_alpha = 1e5,
                              nbsamples_beta = 2e3){

  nbPixValid <- get_valid_pixels_from_tiles(feature_dir, plots, nbCPU = 1,
                                            mask_dir = mask_dir)
  # sample the plot network for alpha diversity
  samples_alpha_terra <- sample_from_plots_alpha(feature_dir = feature_dir,
                                                 list_features = list_features,
                                                 plots = plots,
                                                 nbPixValid = nbPixValid,
                                                 mask_dir = mask_dir,
                                                 nbsamples_alpha = nbsamples_alpha,
                                                 nbCPU = nbCPU)

  # sample the plot network for beta diversity
  samples_beta_terra <- sample_from_plots_beta(feature_dir = feature_dir,
                                               list_features = list_features,
                                               plots = plots,
                                               window_size = window_size,
                                               nbPixValid = nbPixValid,
                                               mask_dir = mask_dir,
                                               nbsamples_beta = nbsamples_beta,
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
