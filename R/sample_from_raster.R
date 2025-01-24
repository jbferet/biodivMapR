#' sample pixels or plots from raster data
#'
#' @param extent_area .
#' @param nbSamples numeric. number of samples to be extracted
#' @param input_rast SpatRaster. raster to extract data from
#' @param input_mask SpatRaster. mask corresponding to raster to extract data from
#' @param window_size numeric. window size for square plots
#' @param capstyle character. shape of the plot, see terra::buffer for info on options for capstyle
#'
#' @return rast_sample dataframe. pixel/plot info extracted from input_rast
#' @importFrom sf st_sample st_as_sf
#' @importFrom terra vect buffer extract res
#' @importFrom stats runif
#' @export

sample_from_raster <- function(extent_area,
                               nbSamples,
                               input_rast,
                               input_mask = NULL,
                               window_size = NULL,
                               capstyle  = 'square'){

  if (nbSamples < 100000 | !is.null(window_size)){
    if (nbSamples > 100000 & !is.null(window_size)){
      message('extracting large number of plots will take some time')
      message('define window_size = NULL to extract pixels only')
    }
    # randomly sample square cells within the image and mask
    samples <- sf::st_sample(x = sf::st_as_sf(extent_area), size = nbSamples)
    # define cell size based on window_size
    if (is.null(window_size)) {
      samples <- terra::vect(samples)
    } else {
      # apply buffer
      bufferSize <- terra::res(input_rast[[1]])[1]*window_size/2
      samples <- terra::buffer(x = terra::vect(samples), width = bufferSize,
                               quadsegs = 8, capstyle  = capstyle)
    }
    rast_sample <- sample_raster(input_rast = input_rast, pix2extract = samples)
    # account for mask if provided
    if (!is.null(input_mask)){
      mask_sample <- terra::extract(x = input_mask, y = samples)
      sel <- which(mask_sample[,2]==1)
      rast_sample <- rast_sample[sel,]
    }
    rast_sample <- clean_NAsInf(rast_sample)
  } else {
    xysamples <- get_xy_samples(input_rast, nbSamples, input_mask = input_mask)
    rast_sample <- sample_raster(input_rast = input_rast,
                                 pix2extract = xysamples,
                                 xy = T)
    rast_sample <- clean_NAsInf(rast_sample)
  }
  return(rast_sample)
}

