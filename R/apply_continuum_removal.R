#' prepares data to run multithreaded continuum removal
#'
#' @param Spectral_Data numeric. initial data matrix (nb samples x nb bands)
#' @param Spectral list. information about spectral bands
#
#' @return samples from image and updated number of pixels to sample if necessary
#' @importFrom snow splitRows
#' @export

apply_continuum_removal <- function(Spectral_Data, Spectral) {
  if (length(Spectral$WaterVapor) > 0) {
    Spectral_Data <- Spectral_Data[, -Spectral$WaterVapor]
  }
  # split data to perform continuum removal on into reasonable amount of data
  nb.Values <- dim(Spectral_Data)[1] * dim(Spectral_Data)[2]
  if (nb.Values > 0) {
    # corresponds to ~ 40 Mb data, but CR tends to requires ~ 10 times memory
    # avoids memory crash
    Max.nb.Values <- 2e6
    nb_CR <- ceiling(nb.Values / Max.nb.Values)
    Spectral_Data <- snow::splitRows(Spectral_Data, nb_CR)
    Spectral_Data_tmp <- lapply(Spectral_Data, FUN = continuumRemoval,
                                Spectral_Bands = Spectral$Wavelength)
    Spectral_Data <- do.call("rbind", Spectral_Data_tmp)
    rm(Spectral_Data_tmp)
  } else {
    # edit 31-jan-2018
    # resize to delete first and last band as in continuum removal
    Spectral_Data <- Spectral_Data[, -c(1, 2)]
  }
  gc()
  return(Spectral_Data)
}
