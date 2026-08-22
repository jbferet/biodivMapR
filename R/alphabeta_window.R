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
#' @return list of alpha and beta diversity metrics corresponding to SSwindow
#' @importFrom stats sd
#' @export

alphabeta_window <- function(SSwindow, nb_clusters, Beta_info, alpha_metrics,
                             beta_metrics, pcelim = 0.02, Hill_order = 1,
                             diss_output = FALSE, p = NULL){
  # get spectral species distribution from individual pixels within a window
  ssd <- lapply(X = SSwindow,FUN = table)
  nb_pix_sunlit <- dim(SSwindow)[1]
  # full spectral species distribution = missing clusters set to 0
  ssd_mat <- do.call('rbind', lapply(X = ssd, FUN = get_normalized_ssd_mat,
                                     nb_clusters = nb_clusters, pcelim = pcelim))

  # get ALPHA diversity
  alpha <- get_alpha(ssd_mat = ssd_mat, alpha_metrics = alpha_metrics,
                     Hill_order = Hill_order)
  alpha <- data.frame(alpha)

  # get BETA diversity
  nb_iter <- nrow(ssd_mat)
  ssd_list <- snow::splitRows(x = ssd_mat, ncl = nb_iter)
  pcoa_diss <- mat_diss <- list()
  if (!is.null(Beta_info) | diss_output){
    # compute dissimilarity with random set
    if (!is.null(Beta_info))
      mat_diss_all <- mapply(FUN = compute_dissimilarity, A = ssd_list,
                             B = Beta_info$ssd,
                             MoreArgs = list(beta_metrics = beta_metrics),
                             SIMPLIFY = FALSE)
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

  list_output <- list('richness_mean' = mean(alpha$Richness),
                      'richness_sd' = stats::sd(alpha$Richness),
                      'shannon_mean' = mean(alpha$Shannon),
                      'shannon_sd' = stats::sd(alpha$Shannon),
                      'simpson_mean' = mean(alpha$Simpson),
                      'simpson_sd' = stats::sd(alpha$Simpson),
                      'hill_mean' = mean(alpha$Hill),
                      'hill_sd' = stats::sd(alpha$Hill),
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
