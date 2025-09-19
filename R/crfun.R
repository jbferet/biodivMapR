#' computes continuum removal for individual spectra
#' @description adapted from propsectr function
#' https://github.com/l-ramirez-lopez/prospectr/blob/main/R/continuumRemoval.R

#' @param x numeric. original spectrum
#' @param wav numeric. central wavelength for the spectral bands
#' @param interpol character. linear or spline
#
#' @return continuum removed spectrum
#' @importFrom grDevices chull
#' @importFrom stats approx splinefun
#' @export

crfun <- function(x, wav, interpol) {
  # need to define close neighbors corresponding to lower and upper bands
  neighbor <- wav[1]/1e4
  # convex hull
  id <- sort(chull(c(wav[1] - neighbor, wav, wav[length(wav)] + neighbor), c(-max(x), x, -max(x))))
  id <- id[-c(1, length(id))] - 1
  cont <- switch(interpol,
                 linear = {
                   approx(x = wav[id], y = x[id], xout = wav, method = "linear")$y
                 },
                 spline = {
                   splinefun(x = wav[id], y = x[id])(wav)
                 }
  )
  return(cont)
}
