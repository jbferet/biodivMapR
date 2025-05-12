#' defines the number of pixels per iteration
#'
#' @param input_rast SpatRaster.
#' @param input_mask SpatRaster.
#' @param nb_pix numeric. maximum number of pixels to extract for kmeans
#' @param nb_iter numeric. nb of iterations averaged to compute diversity indices
#'
#' @return Pix_Per_Iter number of pixels per iteration
#' @importFrom terra values
#' @export

define_pixels_per_iter <- function(input_rast, input_mask = NULL,
                                   nb_iter = 10, nb_pix = 1e5) {
  # check nb of pixels available from image or mask
  nb_pixels_Sunlit <- dim(input_rast)[1] * dim(input_rast)[2]
  if (!is.null(input_mask))
    nb_pixels_Sunlit <- sum(terra::values(input_mask),na.rm = TRUE)
  if (nb_pixels_Sunlit<nb_pix)
    nb_pix <- nb_pixels_Sunlit
  # adjust the number of pixels per iteration
  Pix_Per_Iter <- floor(nb_pix/nb_iter)
  return(Pix_Per_Iter)
}
