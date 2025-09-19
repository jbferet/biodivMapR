#' @title Continuum Removal
#' @description adapted from continuumRemoval function in R package propsectr
#' https://github.com/l-ramirez-lopez/prospectr/blob/main/R/continuumRemoval.R
#'
#' Compute the continuum removed values of a data matrix or vector
#' @usage
#' continuumRemoval(X, wav, type = c("R", "A"),
#'                  interpol = c("linear", "spline"),
#'                  method = c("division", "substraction"))
#' @param X a numeric matrix or vector to process (optionally a data frame that can
#' be coerced to a numerical matrix).
#' @param wav optional. A numeric vector of band positions.
#' @param type the type of data: 'R' for reflectance (default), 'A' for
#' absorbance.
#' @param interpol the interpolation method between points on the convex hull:
#' 'linear' (default) or 'spline'.
#' @param method normalization method: 'division' (default) or 'subtraction'
#' (see details section).
#' @author Antoine Stevens & \href{https://orcid.org/0000-0002-5369-5120}{Leonardo Ramirez-Lopez}
#' @return a matrix or vector with the filtered spectra.
#' @details
#' The continuum removal technique was introduced by Clark and Roush (1984)
#' as a method to highlight energy absorption features of minerals.
#' It can be viewed as a way to perform albedo normalization.
#' The algorithm find points lying on the convex hull (local maxima or envelope)
#' of a spectrum, connects the points by linear or spline interpolation and
#' normalizes the spectrum by dividing (or subtracting) the input data by the
#' interpolated line.
#' @references
#' Clark, R.N., and Roush, T.L., 1984. Reflectance Spectroscopy: Quantitative
#' Analysis Techniques for Remote Sensing Applications. J. Geophys. Res. 89,
#' 6329-6340.
#' @export

continuumRemoval <- function(X,
                             wav,
                             type = c("R", "A"),
                             interpol = c("linear", "spline"),
                             method = c("division", "substraction")) {
  if (is.data.frame(X))
    X <- as.matrix(X)

  type <- match.arg(type)
  interpol <- match.arg(interpol)
  method <- match.arg(method)

  if (type == "A")
    X <- 1 / X

  if (is.matrix(X)) {
    if (missing(wav))
      wav <- seq_len(ncol(X))
    if (length(wav) != ncol(X))
      stop("length(wav) should be equal to ncol(X)")
    cont <- t(apply(X, 1, function(x) crfun(x, wav, interpol)))
  } else {
    cont <- crfun(X, wav, interpol)
  }

  if (method == "division") {
    cr <- X / cont
  } # like ENVI
  else {
    cr <- 1 + X - cont
  }

  if (type == "A")
    cr <- 1 / cr - 1

  if (is.matrix(X)) {
    colnames(cr) <- wav
    rownames(cr) <- rownames(X)
  } else {
    names(cr) <- wav
  }
  return(cr)
}
