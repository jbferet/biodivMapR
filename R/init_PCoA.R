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
#' @import cli
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor handlers with_progress
#' @importFrom parallel makeCluster stopCluster
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

    # # list per plot
    # rast_sample <- rast_sample %>% group_split(ID, .keep = F)
    # if (verbose ==T) message('compute spectral species from beta plots')
    # # compute spectral species for each plot
    # ResDist <- lapply(X = rast_sample, FUN = apply_kmeans,
    #                   Kmeans_info = Kmeans_info,
    #                   SelectBands = SelectBands)
    #
    # # spectral species distribution
    # SSdist <- list()
    # for (iter in names(ResDist[[1]])) SSdist[[iter]] <- lapply(ResDist, '[[',iter)
    # # get nbIter and nbclusters
    # nbIter <- length(Kmeans_info$Centroids)
    # nbclusters <- dim(Kmeans_info$Centroids[[1]])[1]
    # # compute spectral species distribution for each cluster & BC dissimilarity
    # if (verbose ==T) message('compute dissimilarity among plots')
    # # plan(multisession, workers = nbCPU)
    # cl <- parallel::makeCluster(nbCPU)
    # plan("cluster", workers = cl)
    # handlers(global = TRUE)
    # handlers("cli")
    # with_progress({
    #   p <- progressr::progressor(steps = nbIter)
    #   Beta_info <- future.apply::future_lapply(SSdist,
    #                                            FUN = get_BCdiss_from_SSD,
    #                                            nbclusters = nbclusters,
    #                                            pcelim = pcelim, p = p)
    # })
    # parallel::stopCluster(cl)
    # plan(sequential)
    # MatBC_iter <- lapply(Beta_info, '[[','MatBC')
    # SSD <- lapply(Beta_info, '[[','SSD')
    # MatBC <- Reduce('+', MatBC_iter)/nbIter
    # MatBCdist <- stats::as.dist(MatBC, diag = FALSE, upper = FALSE)
    # BetaPCO <- labdsv::pco(MatBCdist, k = dimPCoA)
    # # Beta_Ordination_sel <- BetaPCO$points
    # Beta_info <- list('SSD' = SSD, 'MatBC' = MatBC, 'BetaPCO' = BetaPCO)
    # if (is.null(Beta_info_save)) Beta_info_save <- file.path(output_dir,
    #                                                          'Beta_info.RData')
    # save(Beta_info, file = Beta_info_save)
  }
  return(Beta_info)
}
