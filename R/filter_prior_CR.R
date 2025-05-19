#' Filter data prior to continuum removal:
#' - values are expected to be real reflectance values between 0 and 10000
#' - negative values may occur, so a +100 value is applied to avoid negative
#' - possibly remaining negative values are set to 0
#' - constant spectra are eliminated
#'
#' @param mat_init initial data matrix, n rows = n samples, p cols = p bands
#' @param spectral_bands numeric. central wavelength for the spectral bands
#
#' @return list. updated mat_init
#' @export

filter_prior_cr <- function(mat_init, spectral_bands) {
  # number of samples to be processed
  nb_samples <- nrow(mat_init)
  # make sure there is no negative values
  # mat_init[mat_init<0] <- mat_init + 100.0
  mat_init[mat_init < 0] <- 0
  # eliminate invariant spectra
  null_sd <- apply(mat_init, 1, sd)
  samples_to_keep <- which(!null_sd == 0 & !is.na(null_sd))
  mat_init <- mat_init[samples_to_keep, ]
  # add negative values to the last column and update spectral bands
  mat_init <- cbind(mat_init, -9999)
  nb_bands <- ncol(mat_init)
  spectral_bands <- c(spectral_bands, spectral_bands[nb_bands-1] + 100)
  return(list("mat_init" = mat_init, "spectral_bands" = spectral_bands,
              "nb_samples" = nb_samples, "samples_to_keep" = samples_to_keep))
}
