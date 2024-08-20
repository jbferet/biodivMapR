#' compute alpha and beta diversity metrics from pixel data corresponding to
#' spectral species extracted from a window
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nbclusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param Hill_order numeric. Hill order
#' @param p list. progressor object for progress bar
#'
#' @return list of alpha and beta diversity metrics
#' @importFrom stats sd
#' @export

alphabeta_window <- function(SSwindow, nbclusters,
                             Beta_info, alphametrics,
                             pcelim = 0.02,
                             Hill_order = 1,
                             p = NULL){
  # get spectral species distribution from individual pixels within a window
  SSD <- lapply(X = SSwindow,FUN = table)
  # get ALPHA diversity
  nbPix_Sunlit <- dim(SSwindow)[1]
  alpha <- lapply(X = SSD,
                  FUN = get_alpha_from_SSD,
                  alphametrics = alphametrics,
                  nbPix_Sunlit = nbPix_Sunlit,
                  pcelim = pcelim,
                  Hill_order = Hill_order)
  # get BETA diversity
  # full spectral species distribution = missing clusters set to 0
  SSD_full <- lapply(X = SSD, FUN = get_SSD_full,
                     nbclusters = nbclusters, pcelim = pcelim)
  MatBCtmp <- list()
  nbIter <- length(SSD_full)
  PCoA_BC <- NULL
  if (!is.null(Beta_info)){
    for (i in 1:nbIter) MatBCtmp[[i]] <- list('mat1' = SSD_full[[i]],
                                              'mat2' = Beta_info$SSD[[i]])
    MatBCtmp0 <- lapply(X = MatBCtmp, FUN = compute_BCdiss, pcelim)
    MatBCtmp <- Reduce('+', MatBCtmp0)/nbIter
    PCoA_BC <- compute_NN_from_ordination(MatBC = MatBCtmp, knn = 3,
                                          PCoA_train = Beta_info$BetaPCO$points)
  }
  if (!is.null(p)) p()
  return(list('richness_mean' = mean(unlist(lapply(alpha, '[[', 'richness'))),
              'richness_sd' = stats::sd(unlist(lapply(alpha, '[[', 'richness'))),
              'shannon_mean' = mean(unlist(lapply(alpha, '[[', 'shannon'))),
              'shannon_sd' = stats::sd(unlist(lapply(alpha, '[[', 'shannon'))),
              'simpson_mean' = mean(unlist(lapply(alpha, '[[', 'simpson'))),
              'simpson_sd' = stats::sd(unlist(lapply(alpha, '[[', 'simpson'))),
              'fisher_mean' = mean(unlist(lapply(alpha, '[[', 'fisher'))),
              'fisher_sd' = stats::sd(unlist(lapply(alpha, '[[', 'fisher'))),
              'hill_mean' = mean(unlist(lapply(alpha, '[[', 'hill'))),
              'hill_sd' = stats::sd(unlist(lapply(alpha, '[[', 'hill'))),
              'PCoA_BC' = PCoA_BC))
}
