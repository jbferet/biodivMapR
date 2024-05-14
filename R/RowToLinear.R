#' get max index for each row and convert into linear index
#' @param MM numeric.
#' @param nbi numeric.
#' @param nbj numeric.
#'
#' @return MaxCont
#' @export

RowToLinear <- function(MM, nbi, nbj) {
  adj <- seq_len(nbi)
  MaxCont <- ((MM - 1) * (nbi)) + adj
  return(MaxCont)
}

