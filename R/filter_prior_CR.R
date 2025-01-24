#' Filter data prior to continuum removal:
#' - values are expected to be real reflectance values between 0 and 10000
#' - negative values may occur, so a +100 value is applied to avoid negative
#' - possibly remaining negative values are set to 0
#' - constant spectra are eliminated
#'
#' @param Minit initial data matrix, n rows = n samples, p cols = p spectral bands
#' @param Spectral_Bands numeric. central wavelength for the spectral bands
#
#' @return list. updated Minit
#' @export

filter_prior_CR <- function(Minit, Spectral_Bands) {
  # number of samples to be processed
  nbSamples <- nrow(Minit)
  # make sure there is no negative values
  # Minit[Minit<0] <- Minit + 100.0
  Minit[Minit < 0] <- 0
  # eliminate invariant spectra
  SD <- apply(Minit, 1, sd)
  keep <- which(!SD == 0 & !is.na(SD))
  Minit <- Minit[keep, ]
  nbSamplesUpDate <- nrow(Minit)
  # add negative values to the last column and update spectral bands
  Minit <- cbind(Minit, -9999)
  nbBands <- ncol(Minit)
  Spectral_Bands <- c(Spectral_Bands, Spectral_Bands[nbBands-1] + 100)
  return(list("Minit" = Minit, "Spectral_Bands" = Spectral_Bands,
              "nbSamples" = nbSamples, "SamplesToKeep" = keep))
}
