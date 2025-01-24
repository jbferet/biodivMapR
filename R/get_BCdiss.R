#' Computes BC dissimilarity for a given matrix
#'
#' @param Mat1 numeric. matrix of spectral species distribution
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param p list. progressor object for progress bar
#'
#' @return BCtmp numeric. BC dissimilarity matrix corresponding to Mat1 and Mat2
#' @export

get_BCdiss <- function(Mat1, pcelim = 0.02, p = NULL){
  SSDList <- list(Mat1, Mat1)
  BCtmp <- compute_BCdiss(SSDList, pcelim)
  if (!is.null(p)) p()
  return(BCtmp)
}
