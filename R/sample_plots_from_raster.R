#' sample pixels or plots from raster data
#'
#' @param extent_area .
#' @param nb_samples numeric. number of samples to be extracted
#' @param input_rast SpatRaster. raster to extract data from
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#' @param input_mask SpatRaster. mask corresponding to raster
#' @param window_size numeric. window size for square plots
#'
#' @return rast_sample dataframe. pixel/plot info extracted from input_rast
#' @export

sample_plots_from_raster <- function(extent_area,
                                     nb_samples, input_rast, min_sun = 0.25,
                                     input_mask = NULL, window_size = NULL){
  # sample pixels / plots
  rast_sample <- sample_from_raster(extent_area = extent_area,
                                    nb_samples = nb_samples,
                                    input_rast = input_rast,
                                    window_size = window_size,
                                    input_mask = input_mask)

  # select only data which is not masked
  # eliminate cells with less than XX% valid pixels
  keepBeta <- as.integer(names(which(table(rast_sample$ID)>=min_sun*window_size**2)))
  keepLines <- which(rast_sample$ID %in% keepBeta)
  rast_sample <- rast_sample[keepLines,]
  return(rast_sample)
}

