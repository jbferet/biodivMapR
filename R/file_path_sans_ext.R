#' Re-implementation of \code{\link{file_path_sans_ext}} in \code{tools}. This
#' version can handle "." just before the file extenstion, unlike the original
#' implementation.
#' taken from https://github.com/cbig/zonator
#'
#' @param x Character vector giving file paths.
#' @param compression	Logical: should compression extension '.gz', '.bz2' or
#'   '.xz' be removed first?
#'
#' @return File path without the file extension.
#'
#' @author Joona Lehtomaki \email{joona.lehtomaki@@gmail.com}
#' @export
#'
file_path_sans_ext <- function(x, compression = FALSE) {
  if (compression)
    x <- sub("[.](gz|bz2|xz)$", "", x)
  sub("([^.]+.+)\\.[[:alnum:]]+$", "\\1", x)
}
