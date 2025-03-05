#' define water vapor bands based on spectral sampling of original image
#'
#' @param input_rast character. path of the image
#' @param Excluded_WL numeric. spectral domains corresponding to water vapor absorption
#'
#' @return bands corresponding to atmospheric water absorption domain
#' @importFrom terra names
#' @export

exclude_spectral_domains <- function(input_rast, Excluded_WL = NULL) {
  # definition of water vapor absorption
  if (is.null(Excluded_WL)) {
    Excluded_WL <- data.frame('min' = c(0, 895, 1320, 1780, 2450),
                              'max' = c(400, 1005, 1480, 2040, 2600))
  }
  check_data(input_data = terra::names(input_rast), arguments = list('filter' = 'WL'))
  Bands2Keep <- seq_len(dim(input_rast)[3])
  wl <- WaterVapor <- integer()
  if (!TRUE %in% is.na(as.numeric(terra::names(input_rast)))) {
    wl <- as.numeric(terra::names(input_rast))
    WaterVapor <- c()
    for (w in seq_len(nrow(Excluded_WL))) {
      WaterVapor <- c(WaterVapor, which(wl > Excluded_WL[w, 'min'] & wl < Excluded_WL[w, 'max']))
    }
    if (!length(WaterVapor) == 0) {
      wl <- wl[-WaterVapor]
      Bands2Keep <- Bands2Keep[-WaterVapor]
    }
  }
  return(list("Wavelength" = wl,
              "WaterVapor" = WaterVapor,
              "Bands2Keep" = Bands2Keep))
}
