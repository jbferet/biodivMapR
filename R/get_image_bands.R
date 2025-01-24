#' gets rank of spectral bands in an image
#'
#' @param Spectral_Bands wavelength (nm) of the spectral bands to be found
#' @param wavelength wavelength (nm) of all wavelengths in the image
#'
#' @return rank of all spectral bands of interest in the image and corresponding wavelength
#' and distance to wavelength
#' @export

get_image_bands <- function(Spectral_Bands, wavelength) {
  ImBand <- Distance2WL <- list()
  for (band in names(Spectral_Bands)) {
    ImBand[[band]] <- Distance2WL[[band]] <- NA
    if (!is.na(Spectral_Bands[[band]])){
      Closest_Band <- order(abs(wavelength - Spectral_Bands[[band]]))[1]
      ImBand[[band]] <- Closest_Band
      Distance2WL[[band]] <- abs(wavelength[Closest_Band] - Spectral_Bands[[band]])
    }
  }
  return(list("ImBand" = data.frame(ImBand),
              "Distance2WL" = data.frame(Distance2WL)))
}
