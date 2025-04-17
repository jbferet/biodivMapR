#' initialize PCoA for beta diversity mapping
#'
#' @param input_rast SpatRaster. raster to extract data from
#' @param output_dir character. Path for output directory
#' @param MinSun numeric. minimum percentage of sunlit pixels
#' @param window_size numeric. window size for square plots
#' @param nbSamples numeric. number of samples to be extracted
#' @param Kmeans_info list. obtained from prepare_init_kmeans
#' @param SelectBands numeric. bands selected from input_rast
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param dimPCoA numeric.
#' @param input_mask SpatRaster. mask corresponding to raster to extract data from
#' @param nbCPU numeric. Number of CPUs available
#' @param Beta_info_save character. path where to save Beta_info
#' @param Beta_info_read character. path where to read Beta_info
#' @param verbose boolean. set true for messages
#'
#' @return list including spectral species distribution & BC diss matrix per plot, BetaPCO model
#' @export

init_PCoA <- function(input_rast, output_dir, window_size, Kmeans_info,
                      SelectBands = NULL, input_mask = NULL, nbSamples = 1000,
                      MinSun = 0.25, pcelim = 0.02, dimPCoA = 3, nbCPU = 1,
                      Beta_info_save = NULL, Beta_info_read = NULL, verbose = T){

  # if path for data required to map beta diversity provided
  if (!is.null(Beta_info_read)) {
    if (file.exists(Beta_info_read)) load(Beta_info_read)
    if (!file.exists(Beta_info_read)){
      print_error_message('Beta_info_file_missing')
      Beta_info_read <- NULL
    }
  }
  # if no path for data required to map beta diversity provided
  if (is.null(Beta_info_read)){
    # define raster extent where to randomly sample square plots
    if (verbose ==T) message('sampling plots to prepare for beta diversity mapping')
    extent_area <- get_raster_extent(input_rast[[1]])
    # sample plots for initialization of beta diversity
    rast_sample <- sample_plots_from_raster(extent_area = extent_area,
                                            nbSamples = nbSamples,
                                            input_rast = input_rast,
                                            MinSun = MinSun,
                                            input_mask = input_mask,
                                            window_size = window_size)
    keepBeta <- length(unique(rast_sample$ID))
    # complement sampling based on plots successfully sampled from selection #1
    RateSuccess <- keepBeta/nbSamples
    if (RateSuccess<1){
      # how much should we sample to complement the initial subset?
      # 10% extra plots to minimize risk to undersample
      nbSamplesExtra <- nbSamples-keepBeta
      nb2add <- ceiling(1.1*nbSamplesExtra/RateSuccess)
      rast_sample2 <- sample_plots_from_raster(extent_area = extent_area,
                                               nbSamples = nb2add,
                                               input_rast = input_rast,
                                               MinSun = MinSun,
                                               input_mask = input_mask,
                                               window_size = window_size)
      # eliminate extra samples
      keepBeta <- as.integer(names(which(table(rast_sample2$ID)>=MinSun*window_size**2)))
      keepBeta <- keepBeta[seq_len(min(c(length(keepBeta),nbSamplesExtra)))]
      keepLines <- which(rast_sample2$ID %in% keepBeta)
      # concatenate raster
      rast_sample2 <- rast_sample2[keepLines,]
      rast_sample2$ID <- rast_sample2$ID + nbSamples
      rast_sample <- rbind(rast_sample, rast_sample2)
    }
    Beta_info <- init_PCoA_samples(rast_sample = rast_sample, output_dir = output_dir,
                                   Kmeans_info = Kmeans_info, SelectBands = SelectBands,
                                   pcelim = pcelim, dimPCoA = dimPCoA, nbCPU = nbCPU,
                                   Beta_info_save = Beta_info_save, verbose = verbose)

  }
  return(Beta_info)
}
