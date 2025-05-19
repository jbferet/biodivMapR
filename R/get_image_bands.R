#' gets rank of spectral bands in an image
#'
#' @param spectral_bands wavelength (nm) of the spectral bands to be found
#' @param wavelength wavelength (nm) of all wavelengths in the image
#'
#' @return rank of all spectral bands of interest in the image and corresponding wavelength
#' and distance to wavelength
#' @export

get_image_bands <- function(spectral_bands, wavelength) {
  image_band <- distance_to_wl <- list()
  for (band in names(spectral_bands)) {
    image_band[[band]] <- distance_to_wl[[band]] <- NA
    if (!is.na(spectral_bands[[band]])){
      closest_band <- order(abs(wavelength - spectral_bands[[band]]))[1]
      image_band[[band]] <- closest_band
      distance_to_wl[[band]] <- abs(wavelength[closest_band] - spectral_bands[[band]])
    }
  }
  return(list("image_band" = data.frame(image_band),
              "distance_to_wl" = data.frame(distance_to_wl)))
}
