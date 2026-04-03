#' compute alpha metrics during sfs
#'
#' @param alpha_metrics character
#' @param obs2optimize numeric
#' @param SSValid numeric
#' @param corrMethod character
#' @param nb_clusters numeric
#' @param Hill_order numeric
#' @param pcelim numeric
#' @param nbPlots_total numeric
#'
#' @return assessed value and correlation
#' @export

alpha_metrics_sfs <- function(alpha_metrics, obs2optimize, SSValid, corrMethod,
                              nb_clusters, Hill_order, pcelim, nbPlots_total){

  windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
  windows_per_plot$win_ID <- list(SSValid$win_ID)
  alphabetaIdx_CPU <- lapply(X = windows_per_plot$SSwindow_perCPU,
                             FUN = alphabeta_window_list,
                             nb_clusters = nb_clusters,
                             alpha_metrics = alpha_metrics,
                             Hill_order = Hill_order, pcelim = pcelim)
  alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = FALSE)
  rm(alphabetaIdx_CPU)
  gc()
  # 7- reshape alpha diversity metrics
  IDwindow <- unlist(windows_per_plot$IDwindow_perCPU)
  selcrit <- list('richness'= 'richness_mean', 'shannon' = 'shannon_mean',
                  'simpson' = 'simpson_mean', 'hill' = 'hill_mean')
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

