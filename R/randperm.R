#' performs random permutation for k samples among the vector defined by a
#' original source: package labdsv https://rdrr.io/cran/labdsv/src/R/pco.R
#
#' @param a numeric. vctor of values
#' @param k numeric. number of samples
#
#' @return list
#' @importFrom stats cmdscale
#' @export
#'

randperm <- function(a, k) {
  n <- length(a)
  if (n == 0 || a[1] == 0) return(c())
  if (n == 1) {
    if (floor(a) != ceiling(a) || a < 1)
      stop("Argument 'a' must be a positive integer.")
    n <- a; a <- 1:a
  }
  if (missing(k)) k <- n
  if (k > n)
    stop("'k' must be smaller or equal to 'a' or length of 'a'.")

  m <- sample(a, size = k, replace = FALSE)
  return(m)
}
