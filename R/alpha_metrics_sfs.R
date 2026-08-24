#' compute alpha metrics from list of spectral species distributions 
#' in the context of sequential feature selection
#' 
#' @param alpha_metrics character. names of alpha diversity metrics to compute
#' @param obs2optimize list. reference value for alpha diversity metrics 
#' @param SSValid numeric. 
#' @param corrMethod character. 'pearson' or 'spearman'
#' @param nb_clusters numeric. number of clusters 
#' @param Hill_order numeric. hill order if hill index to be computed
#' @param pcelim numeric. threshold value of relative abundance to keep
#' @param nbPlots_total numeric. total number of plots in corresponding to SSValid
#'
#' @return list. assessed value and correlation with obs2optimize values
#' @export

alpha_metrics_sfs <- function(alpha_metrics, obs2optimize, SSValid, corrMethod,
                              nb_clusters, Hill_order = 1, pcelim, nbPlots_total){

  windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
  windows_per_plot$win_ID <- list(SSValid$win_ID)
  alphabetaIdx <- lapply(X = windows_per_plot$SSwindow_perCPU,
                               # FUN = alphabeta_window_list,
                               FUN = alphabeta_chunk,
                             nb_clusters = nb_clusters,
                             alpha_metrics = alpha_metrics,
                             Hill_order = Hill_order, pcelim = pcelim)
  # alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = FALSE)
  # rm(alphabetaIdx_CPU)
  gc()
  # 7- reshape alpha diversity metrics
  IDwindow <- unlist(windows_per_plot$IDwindow_perCPU)
  selcrit <- list('richness'= 'richness_mean', 'shannon' = 'shannon_mean',
                  'simpson' = 'simpson_mean', 'hill' = 'hill_mean')
  CorrVal <- Assess <- list()
  for (crit in names(selcrit)){
    if (!is.null(obs2optimize[[crit]])){
      Assess[[crit]] <- rep(x = NA, nbPlots_total)
      Assess[[crit]][IDwindow] <- unlist(lapply(alphabetaIdx,'[[',
                                                selcrit[[crit]]))
      CorrVal[[crit]]  <- cor.test(obs2optimize[[crit]],
                                   Assess[[crit]],
                                   method = corrMethod)$estimate
    } else {
      CorrVal[[crit]] <- Assess[[crit]] <- NA
    }
  }
  return(list('Assess' = Assess, 'CorrVal' = CorrVal, 'IDwindow' = IDwindow))
}

