#' get xy of pixels to sample from raster
#'
#' @param input_rast SpatRaster. raster to extract data from
#' @param nbSamples numeric. number of samples to be extracted
#' @param input_mask SpatRaster. mask corresponding to raster to extract data from
#'
#' @return xy matrix
#' @importFrom pracma randperm
#' @importFrom terra values xyFromCell
#' @export

get_xy_samples <- function(input_rast, nbSamples, input_mask = NULL){
  extent_area <- get_raster_extent(input_rast[[1]])
  latlon <- ext(extent_area)
  # check available nb of pixels
  nbPixels <- dim(input_rast[[1]])[1]*dim(input_rast[[1]])[2]
  # adjust if mask provided
  if (!is.null(input_mask)) nbPixels <- sum(terra::values(input_mask),na.rm = T)
  if (nbSamples>nbPixels) nbSamples <- nbPixels
  choicePix <- pracma::randperm(seq_len(nbPixels),nbSamples)
  if (!is.null(input_mask)){
    whichPix <- which(terra::values(input_mask)==1)
  } else {
    whichPix <- seq_len(nbPixels)
  }
  choicePix <- whichPix[choicePix]
  xycol <- terra::xyFromCell(object = input_rast[[1]],choicePix)
  # lon <- stats::runif(n = nbSamples,min = latlon[1L], max = latlon[2L])
  # lat <- stats::runif(n = nbSamples,min = latlon[3L], max = latlon[4L])
  # latlon <- data.frame('lon' = lon, 'lat' = lat)
  return(xycol)
}
