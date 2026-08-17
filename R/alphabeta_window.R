#' compute alpha and beta diversity metrics from pixel data corresponding to
#' spectral species extracted from a window
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param beta_metrics list. beta diversity metrics
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#' @param Hill_order numeric. Hill order
#' @param diss_output boolean set TRUE to get dissimilarity matrices as outputs
#' @param p list. progressor object for progress bar
#'
#' @return list of alpha and beta diversity metrics
#' @importFrom stats sd
#' @export

alphabeta_window <- function(SSwindow, nb_clusters, Beta_info, alpha_metrics,
                             beta_metrics, pcelim = 0.02, Hill_order = 1,
                             diss_output = FALSE, p = NULL){
  # get spectral species distribution from individual pixels within a window
  ssd <- lapply(X = SSwindow,FUN = table)
  # get ALPHA diversity
  nb_pix_sunlit <- dim(SSwindow)[1]
  alpha <- lapply(X = ssd,
                  FUN = get_alpha_from_ssd,
                  alpha_metrics = alpha_metrics,
                  nb_pix_sunlit = nb_pix_sunlit,
                  pcelim = pcelim,
                  hill_order = Hill_order)
  # get BETA diversity
  # full spectral species distribution = missing clusters set to 0
  ssd_full <- lapply(X = ssd, FUN = get_normalized_ssd,
                     nb_clusters = nb_clusters, pcelim = pcelim)
  nb_iter <- length(ssd_full)
  pcoa_diss <- mat_diss <- list()

  if (!is.null(Beta_info) | diss_output){
    # for (beta in beta_metrics){
    # for (i in seq_len(nb_iter))
    #   mat_bc[[i]] <- list('mat1' = ssd_full[[i]],
    #                       'mat2' = Beta_info$SSD[[i]])
    # compute dissimilarity with random set 
    if (!is.null(Beta_info)){
      mat_diss_all <- mapply(FUN = compute_dissimilarity, A = ssd_full,
                             B = Beta_info$ssd,
                             MoreArgs = list(beta_metrics = beta_metrics),
                             SIMPLIFY = FALSE)
    }
    # # compute dissimilarity with self when testing optimal nb of clusters 
    # if (diss_output){
    #   mat_diss_all <- mapply(FUN = compute_dissimilarity, A = ssd_full,
    #                          B = ssd_full,
    #                          MoreArgs = list(beta_metrics = beta_metrics),
    #                          SIMPLIFY = FALSE)
    # }
    
    for (beta in beta_metrics){
      mat_diss_ind <- lapply(X = mat_diss_all, FUN = '[[', beta)
      mean_mat_diss <- Reduce('+', mat_diss_ind)/nb_iter
      if (!is.null(Beta_info$beta_pco)){
        pcoa_diss[[beta]] <- compute_nn_from_ordination(
          mat_bc = mean_mat_diss,
          knn = 3,
          pcoa_train = Beta_info$beta_pco[[beta]]$points)
      }
      if (diss_output){
        mat_diss[[beta]] <- mean_mat_diss
      }
    }
  }
  if (!is.null(p))
    p()

  list_output <- list('richness_mean' = mean(unlist(lapply(alpha, '[[', 'richness'))),
                      'richness_sd' = stats::sd(unlist(lapply(alpha,'[[','richness'))),
                      'shannon_mean' = mean(unlist(lapply(alpha, '[[', 'shannon'))),
                      'shannon_sd' = stats::sd(unlist(lapply(alpha, '[[', 'shannon'))),
                      'simpson_mean' = mean(unlist(lapply(alpha, '[[', 'simpson'))),
                      'simpson_sd' = stats::sd(unlist(lapply(alpha, '[[', 'simpson'))),
                      'fisher_mean' = mean(unlist(lapply(alpha, '[[', 'fisher'))),
                      'fisher_sd' = stats::sd(unlist(lapply(alpha, '[[', 'fisher'))),
                      'hill_mean' = mean(unlist(lapply(alpha, '[[', 'hill'))),
                      'hill_sd' = stats::sd(unlist(lapply(alpha, '[[', 'hill'))),
                      'pcoa_bray' = pcoa_diss[['bray']],
                      'pcoa_brayturn' = pcoa_diss[['brayturn']],
                      'pcoa_simpson_diss' = pcoa_diss[['simpson_diss']],
                      'pcoa_jaccard' = pcoa_diss[['jaccard']],
                      'pcoa_jaccardturn' = pcoa_diss[['jaccardturn']],
                      'pcoa_sorensen' = pcoa_diss[['sorensen']])
  if (diss_output)
    for (beta in beta_metrics)
      list_output[[paste0('matrix_', beta)]] <- mat_diss[[beta]]
  return(list_output)
}
