#' explore disimilarity metrics for a range of number of clusters
#'
#' @param input_rast SpatRaster. raster to extract data from
#' @param output_dir character. Path for output directory
#' @param min_sun numeric. minimum percentage of sunlit pixels
#' @param window_size numeric. window size for square plots
#' @param nb_samples numeric. number of samples to be extracted
#' @param nbClust_list numeric. number of clusters to be tested
#' @param Kmeans_info list. obtained from prepare_init_kmeans
#' @param selected_bands numeric. bands selected from input_rast
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param dimPCoA numeric.
#' @param input_mask SpatRaster. mask corresponding to raster to extract data from
#' @param nbCPU numeric. Number of CPUs available
#' @param beta_metrics character. name of beta dissimilarity metrics
#' Available options are "bray", "brayturn", "jaccard", "jaccardturn",
#' "sorensen", "simpson_diss"
#' @param verbose boolean. set true for messages
#'
#' @return list including spectral species distribution & BC diss matrix per plot, BetaPCO model
#' @export

explore_dissimilarity <- function(
    input_rast, output_dir, window_size, Kmeans_info, selected_bands = NULL, 
    input_mask = NULL, nb_samples = 1000, nbClust_list = 50, min_sun = 0.25, 
    pcelim = 0.02, dimPCoA = 3, nbCPU = 1, beta_metrics = 'bray', verbose = TRUE){
  
  # if no path for data required to map beta diversity provided
  # define raster extent where to randomly sample square plots
  if (verbose)
    message('sampling plots to prepare for beta diversity mapping')
  extent_area <- get_raster_extent(input_rast[[1]])
  # sample plots for initialization of beta diversity
  rast_sample <- sample_plots_from_raster(extent_area = extent_area,
                                          nb_samples = nb_samples,
                                          input_rast = input_rast,
                                          min_sun = min_sun,
                                          input_mask = input_mask,
                                          window_size = window_size)
  keepBeta <- length(unique(rast_sample$ID))
  # complement sampling based on plots successfully sampled from selection #1
  RateSuccess <- keepBeta/nb_samples
  if (RateSuccess<1){
    # how much should we sample to complement the initial subset?
    # 10% extra plots to minimize risk to undersample
    nb_samplesExtra <- nb_samples-keepBeta
    nb2add <- ceiling(1.1*nb_samplesExtra/RateSuccess)
    rast_sample2 <- sample_plots_from_raster(extent_area = extent_area,
                                             nb_samples = nb2add,
                                             input_rast = input_rast,
                                             min_sun = min_sun,
                                             input_mask = input_mask,
                                             window_size = window_size)
    # eliminate extra samples
    keepBeta <- as.integer(names(which(table(rast_sample2$ID)>=min_sun*window_size**2)))
    keepBeta <- keepBeta[seq_len(min(c(length(keepBeta),nb_samplesExtra)))]
    keepLines <- which(rast_sample2$ID %in% keepBeta)
    # concatenate raster
    rast_sample2 <- rast_sample2[keepLines,]
    rast_sample2$ID <- rast_sample2$ID + nb_samples
    rast_sample <- rbind(rast_sample, rast_sample2)
  }
  
  
  
  if (nbCPU>1){
    registerDoFuture()
    cl <- parallel::makeCluster(nbCPU)
    with(plan("cluster", workers = cl), local = TRUE)
    nbclust <- NULL
    get_kmeans_list <- function() {
      foreach(nbclust = seq_along(nbClust_list)) %dopar% {
        Beta_info_save <- file.path(output_dir, paste0('Beta_info_',nbClust_list[nbclust],'.RData'));
        Beta_info <- init_PCoA_samples(rast_sample = rast_sample,
                                       output_dir = output_dir,
                                       Kmeans_info = Kmeans_info[[nbclust]],
                                       selected_bands = selected_bands,
                                       pcelim = pcelim, dimPCoA = dimPCoA,
                                       nbCPU = nbCPU,
                                       beta_metrics = beta_metrics,
                                       Beta_info_save = Beta_info_save,
                                       verbose = verbose)
        return(Beta_info)
      }
    }
    kmeans_info <- get_kmeans_list()
    plan(sequential)
  } else {
    Beta_info <- lapply(X = seq_along(nbClust_list),
                        function(x) {
                          Beta_info_save <- file.path(output_dir, paste0('Beta_info_',nbClust_list[x],'.RData'));
                          init_PCoA_samples(rast_sample = rast_sample,
                                            output_dir = output_dir,
                                            Kmeans_info = Kmeans_info[[x]],
                                            selected_bands = selected_bands,
                                            pcelim = pcelim, dimPCoA = dimPCoA,
                                            nbCPU = nbCPU,
                                            beta_metrics = beta_metrics,
                                            Beta_info_save = Beta_info_save,
                                            verbose = FALSE)
                        })
  }
  return(Beta_info)
}
