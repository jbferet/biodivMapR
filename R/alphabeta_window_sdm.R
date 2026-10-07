#' compute alpha and beta diversity metrics from pixel data corresponding to
#' spectral species extracted from a window
#'
#' @param ssd list.
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param beta_metrics list. beta diversity metrics
#' @param Hill_order numeric. Hill order
#' @param nbCPU numeric. Number of CPUs available
#'
#' @return list of alpha and beta diversity metrics corresponding to ssd
#' @import cli
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom parallel makeCluster stopCluster
#' @importFrom dplyr filter %>%
#' @export

alphabeta_window_sdm <- function(ssd, Beta_info, alpha_metrics, beta_metrics,
                                 Hill_order = 1, nbCPU = 1){
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


  get_pcoa_from_bc <- function(ssd, Beta_info, beta_metrics, p = NULL){

    mat_diss_all <- compute_dissimilarity(A = matrix(ssd, nrow = 1),
                                          B = Beta_info$ssd,
                                          beta_metrics = beta_metrics)


    # mat_bc_tmp <- dissUtils::diss(matrix(ssd, nrow = 1), Beta_info$SSD,
    #                               method = 'braycurtis')
    # mat_bc <- list('mat1' = matrix(ssd, nrow = 1),
    #                'mat2' = Beta_info$SSD)
    # mat_bc_tmp <- compute_bc_diss(ssd_list = mat_bc, pcelim = 0)
    pcoa_diss <- list()
    for (beta in beta_metrics)
      pcoa_diss[[beta]] <- compute_nn_from_ordination(mat_bc = mat_diss_all[[beta]], knn = 3,
                                            pcoa_train = Beta_info$beta_pco[[beta]]$points)
    if (!is.null(p))
      p()
    return(pcoa_diss)
  }

  message('compute beta diversity')
  if (nbCPU == 1){
    pcoa_diss <- if (requireNamespace("pbapply", quietly = TRUE)) {
      pbapply::pblapply(X = ssd,
                        FUN = get_pcoa_from_bc,
                        Beta_info = Beta_info,
                        beta_metrics = beta_metrics)
    } else {
      lapply(X = ssd,
             FUN = get_pcoa_from_bc,
             Beta_info = Beta_info,
             beta_metrics = beta_metrics)
    }
  } else {
    cl <- parallel::makeCluster(nbCPU)
    parallel::clusterEvalQ(cl, {library(biodivMapR)})
    with(plan("cluster", workers = cl), local = TRUE)
    funct <- future.apply::future_lapply
    future.seed = TRUE
    pcoa_diss <- funct(X = ssd,
                       FUN = get_pcoa_from_bc,
                       Beta_info = Beta_info,
                       beta_metrics = beta_metrics,
                       future.seed = future.seed,
                       future.chunk.size = NULL,
                       future.scheduling = 1)
    parallel::stopCluster(cl)
    plan(sequential)
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
              'pcoa_bray' = do.call(rbind,lapply(pcoa_diss, '[[', 'bray')),
              'pcoa_brayturn' = do.call(rbind,lapply(pcoa_diss, '[[', 'brayturn')),
              'pcoa_simpson_diss' = do.call(rbind,lapply(pcoa_diss, '[[', 'simpson_diss')),
              'pcoa_jaccard' = do.call(rbind,lapply(pcoa_diss, '[[', 'jaccard')),
              'pcoa_jaccardturn' = do.call(rbind,lapply(pcoa_diss, '[[', 'jaccardturn')),
              'pcoa_sorensen' = do.call(rbind,lapply(pcoa_diss, '[[', 'sorensen'))))
}
