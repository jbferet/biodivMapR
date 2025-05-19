#' compute alpha and beta diversity metrics from pixel data corresponding to
#' spectral species extracted from a window
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alphametrics list. alpha diversity metrics
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#' @param Hill_order numeric. Hill order
#' @param p list. progressor object for progress bar
#'
#' @return list of alpha and beta diversity metrics
#' @importFrom stats sd
#' @export

alphabeta_window <- function(SSwindow, nb_clusters,
                             Beta_info, alphametrics,
                             pcelim = 0.02,
                             Hill_order = 1,
                             p = NULL){
  # get spectral species distribution from individual pixels within a window
  ssd <- lapply(X = SSwindow,FUN = table)
  # get ALPHA diversity
  nb_pix_sunlit <- dim(SSwindow)[1]
  alpha <- lapply(X = ssd,
                  FUN = get_alpha_from_ssd,
                  alphametrics = alphametrics,
                  nb_pix_sunlit = nb_pix_sunlit,
                  pcelim = pcelim,
                  hill_order = Hill_order)
  # get BETA diversity
  # full spectral species distribution = missing clusters set to 0
  ssd_full <- lapply(X = ssd, FUN = get_ssd_full,
                     nb_clusters = nb_clusters, pcelim = pcelim)
  mat_bc <- list()
  nb_iter <- length(ssd_full)
  pcoa_bc <- NULL
  if (!is.null(Beta_info)){
    for (i in seq_len(nb_iter))
      mat_bc[[i]] <- list('mat1' = ssd_full[[i]],
                          'mat2' = Beta_info$SSD[[i]])
    mat_bc_tmp <- lapply(X = mat_bc, FUN = compute_bc_diss, pcelim)
    mat_bc <- Reduce('+', mat_bc_tmp)/nb_iter
    pcoa_bc <- compute_nn_from_ordination(mat_bc = mat_bc, knn = 3,
                                          pcoa_train = Beta_info$BetaPCO$points)
  }
  if (!is.null(p))
    p()
  return(list('richness_mean' = mean(unlist(lapply(alpha, '[[', 'richness'))),
              'richness_sd' = stats::sd(unlist(lapply(alpha,'[[','richness'))),
              'shannon_mean' = mean(unlist(lapply(alpha, '[[', 'shannon'))),
              'shannon_sd' = stats::sd(unlist(lapply(alpha, '[[', 'shannon'))),
              'simpson_mean' = mean(unlist(lapply(alpha, '[[', 'simpson'))),
              'simpson_sd' = stats::sd(unlist(lapply(alpha, '[[', 'simpson'))),
              'fisher_mean' = mean(unlist(lapply(alpha, '[[', 'fisher'))),
              'fisher_sd' = stats::sd(unlist(lapply(alpha, '[[', 'fisher'))),
              'hill_mean' = mean(unlist(lapply(alpha, '[[', 'hill'))),
              'hill_sd' = stats::sd(unlist(lapply(alpha, '[[', 'hill'))),
              'PCoA_BC' = pcoa_bc))
}
