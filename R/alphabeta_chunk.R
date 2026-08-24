#' apply alphabeta_window to a list of lists
#'
#' @param SSwindow dataframe. spectral species corresponding to a raster subset
#' (window = elementary spatial unit of process)
#' @param nb_clusters numeric. number of clusters used in kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param alpha_metrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#' @param beta_metrics list. beta diversity metrics
#' @param diss_output boolean set TRUE to get dissimilarity matrices as outputs
#' @param pcelim numeric. min proportion of pixels to consider spectral species
#' @param p list. progressor object for progress bar
#' @param ... list. additional parameters
#'
#' @return list of alpha and beta diversity metrics
#'
#' @importFrom dplyr summarize group_by across mutate where %>%
#' @importFrom tidyr nest
#' @importFrom purrr map
#' @export

alphabeta_chunk <- function(
    SSwindow, nb_clusters, Beta_info = NULL, alpha_metrics = 'shannon',
    Hill_order = 1, beta_metrics = 'bray', diss_output = FALSE, pcelim = 0.02,
    p = NULL, ...){


  # ssd <- lapply(X = lapply(X = SSwindow,FUN = table),FUN = table)
  nb_windows <- length(SSwindow)
  nb_pix_sunlit <- dim(SSwindow[[1]])[1]
  # full spectral species distribution = missing clusters set to 0
  ssd_mat <- mapply(ssd_list = SSwindow, windows_id = seq_len(nb_windows),
                    FUN = get_normalized_ssd_list,
                    MoreArgs = list(nb_clusters = nb_clusters, pcelim = pcelim),
                    SIMPLIFY = FALSE)
  ssd_mat <- as.data.frame(do.call('rbind', ssd_mat))
  window_id <- ssd_mat$windows_id
  iter <- ssd_mat$iter
  nb_iter <- length(unique(iter))
  ssd_mat$windows_id <- ssd_mat$iter <- NULL
  ssd_mat_alpha <- as.matrix(ssd_mat)

  # compute alpha diversity
  alpha <- get_alpha(ssd_mat = ssd_mat_alpha, alpha_metrics = alpha_metrics,
                     Hill_order = Hill_order)
  alpha <- data.frame(alpha)
  alpha$ID <- window_id

  alpha_mean <- alpha %>%
    group_by(ID) %>%
    summarize(across(where(is.numeric), function(x) mean(x, na.rm = TRUE)))

  alpha_sd <- alpha %>%
    group_by(ID) %>%
    summarize(across(where(is.numeric), function(x) sd(x, na.rm = TRUE)))

  # compute beta diversity
  pcoa_diss <- mat_diss <- list()
  if (!is.null(Beta_info) | diss_output){
    # group normalized ssd by iter & convert nested dataframe to matrix
    ssd_mat$iter <- iter
    ssd_beta <- ssd_mat %>%
      group_by(iter)  %>%
      nest() %>%
      mutate(
        data_matrix = map(data, ~ as.matrix(.x))
      )
    # compute dissimilarity with random set
    if (!is.null(Beta_info))
      mat_diss_all <- mapply(FUN = compute_dissimilarity, A = ssd_beta$data_matrix,
                             B = Beta_info$ssd,
                             MoreArgs = list(beta_metrics = beta_metrics),
                             SIMPLIFY = FALSE)
    for (beta in beta_metrics){
      mat_diss_ind <- lapply(X = mat_diss_all, FUN = '[[', beta)
      mean_mat_diss <- Reduce('+', mat_diss_ind)/nb_iter
      if (!is.null(Beta_info$beta_pco))
        pcoa_diss[[beta]] <- compute_nn_from_ordination(
          mat_bc = mean_mat_diss,
          knn = 3,
          pcoa_train = Beta_info$beta_pco[[beta]]$points)
      if (diss_output)
        mat_diss[[beta]] <- mean_mat_diss
    }
  }
  if (!is.null(p))
    p()


  list_output <- list('richness_mean' = alpha_mean$Richness,
                      'richness_sd' = alpha_sd$Richness,
                      'shannon_mean' = alpha_mean$Shannon,
                      'shannon_sd' = alpha_sd$Shannon,
                      'simpson_mean' = alpha_mean$Simpson,
                      'simpson_sd' = alpha_sd$Simpson,
                      'hill_mean' = alpha_mean$Hill,
                      'hill_sd' = alpha_sd$Hill,
                      'pcoa_bray' = pcoa_diss[['bray']],
                      'pcoa_brayturn' = pcoa_diss[['brayturn']],
                      'pcoa_simpson_diss' = pcoa_diss[['simpson_diss']],
                      'pcoa_jaccard' = pcoa_diss[['jaccard']],
                      'pcoa_jaccardturn' = pcoa_diss[['jaccardturn']],
                      'pcoa_sorensen' = pcoa_diss[['sorensen']])
  if (diss_output)
    for (beta in beta_metrics)
      list_output[[paste0('matrix_', beta)]] <- mat_diss[[beta]]

  # alphabeta_idx <- lapply(X = SSwindow,
  #                         FUN = alphabeta_window,
  #                         nb_clusters = nb_clusters,
  #                         Beta_info = Beta_info,
  #                         alpha_metrics = alpha_metrics,
  #                         Hill_order = Hill_order,
  #                         beta_metrics = beta_metrics,
  #                         diss_output = diss_output,
  #                         pcelim = pcelim)
  return(list_output)
}
