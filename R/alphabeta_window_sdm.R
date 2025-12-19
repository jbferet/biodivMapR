#' compute alpha and beta diversity metrics from pixel data corresponding to
#' spectral species extracted from a window
#'
#' @param ssd list.
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#'
#' @return list of alpha and beta diversity metrics
#' @importFrom pbapply pblapply
#' @importFrom dissUtils diss
#' @export

alphabeta_window_sdm <- function(ssd, Beta_info, alpha_metrics, Hill_order = 1){
  # get ALPHA diversity
  nb_pix_sunlit <- length(ssd[[1]])
  message('compute alpha diversity')
  alpha <- if (requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply(X = ssd,
                      FUN = get_alpha_from_ssd,
                      alpha_metrics = alpha_metrics,
                      nb_pix_sunlit = nb_pix_sunlit,
                      pcelim = 0,
                      hill_order = Hill_order)
  } else {
    lapply(X = ssd,
           FUN = get_alpha_from_ssd,
           alpha_metrics = alpha_metrics,
           nb_pix_sunlit = nb_pix_sunlit,
           pcelim = 0,
           hill_order = Hill_order)
  }

  # get BETA diversity
  # full spectral species distribution = missing clusters set to 0
  for (i in 1:length(ssd)){
    if (any(is.na(ssd[[i]])))
      ssd[[i]][which(is.na(ssd[[i]]))] <- 0
    ssd[[i]] <- ssd[[i]]/sum(ssd[[i]])
  }


  get_poca_from_bc <- function(ssd, Beta_info, p = NULL){
    mat_bc_tmp <- dissUtils::diss(matrix(ssd, nrow = 1), Beta_info$SSD,
                                  method = 'braycurtis')
    # mat_bc <- list('mat1' = matrix(ssd, nrow = 1),
    #                'mat2' = Beta_info$SSD)
    # mat_bc_tmp <- compute_bc_diss(ssd_list = mat_bc, pcelim = 0)
    pcoa_bc <- compute_nn_from_ordination(mat_bc = mat_bc_tmp, knn = 3,
                                          pcoa_train = Beta_info$BetaPCO$points)
    if (!is.null(p))
      p()
    return(pcoa_bc)
  }

  message('compute beta diversity')
  pcoa_bc <- if (requireNamespace("pbapply", quietly = TRUE)) {
    pbapply::pblapply(X = ssd,
                      FUN = get_poca_from_bc,
                      Beta_info = Beta_info)
  } else {
    lapply(X = ssd,
           FUN = get_poca_from_bc,
           Beta_info = Beta_info)
  }
  # for (i in seq_along(ssd)){
  #   mat_bc <- list('mat1' = matrix(ssd[[i]], nrow = 1),
  #                  'mat2' = Beta_info$SSD)
  #   mat_bc_tmp <- compute_bc_diss(ssd_list = mat_bc, pcelim = pcelim)
  #   pcoa_bc[[i]] <- compute_nn_from_ordination(mat_bc = mat_bc_tmp, knn = 3,
  #                                              pcoa_train = Beta_info$BetaPCO$points)
  # }
  return(list('richness' = unlist(lapply(alpha, '[[', 'richness')),
              'shannon' = unlist(lapply(alpha, '[[', 'shannon')),
              'simpson' = unlist(lapply(alpha, '[[', 'simpson')),
              'hill' = unlist(lapply(alpha, '[[', 'hill')),
              'PCoA_BC' = pcoa_bc))
}
