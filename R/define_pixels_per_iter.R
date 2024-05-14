#' defines the number of pixels per iteration
#'
#' @param input_rast SpatRaster.
#' @param input_mask SpatRaster.
#' @param nbPix numeric. maximum number of pixels to extract for kmeans
#' @param nbIter numeric. nb of iterations averaged to compute diversity indices
#'
#' @return Pix_Per_Iter number of pixels per iteration
#' @importFrom terra values
#' @export

define_pixels_per_iter <- function(input_rast, input_mask = NULL,
                                   nbIter = 20, nbPix = 100000) {
  # check nb of pixels available from image or mask
  nbPixels_Sunlit <- dim(input_rast)[1] * dim(input_rast)[2]
  if (!is.null(input_mask)) nbPixels_Sunlit <- sum(terra::values(input_mask),na.rm = T)
  if (nbPixels_Sunlit<nbPix) nbPix <- nbPixels_Sunlit
  # adjust the number of pixels per iteration
  Pix_Per_Iter <- floor(nbPix/nbIter)
  return(Pix_Per_Iter)
}
