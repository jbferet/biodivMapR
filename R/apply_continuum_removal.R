#' prepares data into list to apply continuum removal
#'
#' @param spectral_data numeric. initial data matrix (nb samples x nb bands)
#' @param spectral list. information about spectral bands
#
#' @return samples from image and updated number of pixels to sample
#' @importFrom snow splitRows
#' @export

apply_continuum_removal <- function(spectral_data, spectral) {
  if (length(spectral$WaterVapor) > 0)
    spectral_data <- spectral_data[, -spectral$WaterVapor]

  # split data to perform continuum removal on into reasonable amount of data
  nb_values <- dim(spectral_data)[1] * dim(spectral_data)[2]
  if (nb_values > 0) {
    # corresponds to ~ 40 Mb data, but CR tends to requires ~ 10 times memory
    # avoids memory crash
    max_nb_values <- 2e6
    nb_cr <- ceiling(nb_values / max_nb_values)
    spectral_data <- snow::splitRows(spectral_data, nb_cr)
    spectral_data_tmp <- lapply(spectral_data, FUN = continuum_removal,
                                spectral_bands = spectral$Wavelength)
    spectral_data <- do.call("rbind", spectral_data_tmp)
    rm(spectral_data_tmp)
  } else {
    # edit 31-jan-2018
    # resize to delete first and last band as in continuum removal
    spectral_data <- spectral_data[, -c(1, 2)]
  }
  gc()
  return(spectral_data)
}
