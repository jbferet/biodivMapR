#' compute alpha and beta diversity metrics from pixel data corresponding to
#' spectral species extracted from a window
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#' @param Hill_order numeric. Hill order
#'
#' @return list of alpha and beta diversity metrics
#' @importFrom dissUtils diss
#' @export

alphabeta_window_classif <- function(SSwindow, nb_clusters,
                                     Beta_info = NULL, alpha_metrics,
                                     pcelim = 0.02,
                                     Hill_order = 1){
  # get spectral species distribution from individual pixels within a window
  ssd <- lapply(X = SSwindow,FUN = table)
  # get ALPHA diversity
  nb_pix_sunlit <- length(SSwindow[[1]])
  alpha <- lapply(X = ssd,
                  FUN = get_alpha_from_ssd,
                  alpha_metrics = alpha_metrics,
                  nb_pix_sunlit = nb_pix_sunlit,
                  pcelim = pcelim,
                  hill_order = Hill_order)

  pcoa_bc <- NULL
  if (!is.null(Beta_info)){
    # get BETA diversity
    # full spectral species distribution = missing clusters set to 0
    ssd_full <- lapply(X = ssd, FUN = get_normalized_ssd,
                       nb_clusters = nb_clusters, pcelim = pcelim)
    mat_bc <- list()
    pcoa_bc <- list()
    for (i in seq_along(ssd_full)){
      mat_bc_tmp <- dissUtils::diss(ssd_full[[i]], Beta_info$SSD,
                                    method = 'braycurtis')
      # mat_bc <- list('mat1' = ssd_full[[i]],
      #                'mat2' = Beta_info$SSD)
      # mat_bc_tmp <- compute_bc_diss(ssd_list = mat_bc, pcelim = pcelim)
      pcoa_bc[[i]] <- compute_nn_from_ordination(mat_bc = mat_bc_tmp, knn = 3,
                                                 pcoa_train = Beta_info$BetaPCO$points)
    }
  }
  return(list('richness' = unlist(lapply(alpha, '[[', 'richness')),
              'shannon' = unlist(lapply(alpha, '[[', 'shannon')),
              'simpson' = unlist(lapply(alpha, '[[', 'simpson')),
              'hill' = unlist(lapply(alpha, '[[', 'hill')),
              'PCoA_BC' = pcoa_bc))
}
