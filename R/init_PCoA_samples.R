#' initialize PCoA for beta diversity mapping based on samples extracted from images
#'
#' @param rast_sample data frame containing samples to use
#' @param output_dir character. Path for output directory
#' @param Kmeans_info list. obtained from prepare_init_kmeans
#' @param SelectBands numeric. bands selected from input_rast
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param dimPCoA numeric.
#' @param nbCPU numeric. Number of CPUs available
#' @param Beta_info_save character. path where to save Beta_info
#' @param verbose boolean. set true for messages
#'
#' @return list including spectral species distribution & BC diss matrix per plot, BetaPCO model
#' @import cli
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor handlers with_progress
#' @importFrom dplyr group_split
#' @importFrom stats as.dist
#' @importFrom parallel makeCluster stopCluster
#' @export

init_PCoA_samples <- function(rast_sample, output_dir, Kmeans_info,
                              SelectBands = NULL, pcelim = 0.02, dimPCoA = 3,
                              nbCPU = 1, Beta_info_save = NULL, verbose = T){

  rast_sample <- clean_NAsInf(rast_sample)
  ID <- NULL
  # list per plot
  rast_sample <- rast_sample %>% group_split(ID, .keep = F)
  if (verbose ==T) message('compute spectral species from beta plots')
  # compute spectral species for each plot
  ResDist <- lapply(X = rast_sample, FUN = apply_kmeans,
                    Kmeans_info = Kmeans_info,
                    SelectBands = SelectBands)

  # spectral species distribution
  SSdist <- list()
  for (iter in names(ResDist[[1]])) SSdist[[iter]] <- lapply(ResDist, '[[',iter)
  # get nbIter and nbclusters
  nbIter <- length(Kmeans_info$Centroids)
  nbclusters <- dim(Kmeans_info$Centroids[[1]])[1]
  # compute spectral species distribution for each cluster & BC dissimilarity
  if (verbose ==T) message('compute dissimilarity among plots')
  # plan(multisession, workers = nbCPU)
  if (nbCPU>1){
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)
    handlers(global = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nbIter)
      Beta_info <- future.apply::future_lapply(SSdist,
                                               FUN = get_BCdiss_from_SSD,
                                               nbclusters = nbclusters,
                                               pcelim = pcelim, p = p,
                                               future.seed = TRUE)
    })
    parallel::stopCluster(cl)
    plan(sequential)
  } else {
    Beta_info <- lapply(X = SSdist,
                        FUN = get_BCdiss_from_SSD,
                        nbclusters = nbclusters,
                        pcelim = pcelim)
  }
  MatBC_iter <- lapply(Beta_info, '[[','MatBC')
  SSD <- lapply(Beta_info, '[[','SSD')
  MatBC <- Reduce('+', MatBC_iter)/nbIter
  MatBCdist <- stats::as.dist(MatBC, diag = FALSE, upper = FALSE)
  BetaPCO <- pco(MatBCdist, k = dimPCoA)
  # Beta_Ordination_sel <- BetaPCO$points
  Beta_info <- list('SSD' = SSD, 'MatBC' = MatBC, 'BetaPCO' = BetaPCO)
  if (is.null(Beta_info_save)) Beta_info_save <- file.path(output_dir,
                                                           'Beta_info.RData')
  save(Beta_info, file = Beta_info_save)
  return(Beta_info)
}
