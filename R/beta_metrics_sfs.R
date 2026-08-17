#' compute beta metrics during sfs
#'
#' @param obs2optimize numeric. number of clusters used in kmeans
#' @param SSValid numeric
#' @param beta_metrics character. list of beta metrics to compute
#' @param nb_clusters numeric
#' @param pcelim numeric
#' @param nb_iter numeric
#' @param nbPlots_total numeric
#' @param Attributes list
#' @param corrMethod character
#' @param IDwindow numeric
#'
#' @return assessed value and correlation
#' @export

beta_metrics_sfs <- function(obs2optimize, SSValid, beta_metrics, nb_clusters, pcelim, nb_iter,
                     nbPlots_total, Attributes, corrMethod, IDwindow){
  win_ID <- NULL
  if (length(beta_metrics)>0){
    SSValid_win <- SSValid %>% group_split(win_ID, .keep = FALSE)
    # spectral species distribution
    SSdist <- list()
    for (iter in names(SSValid_win[[1]]))
      SSdist[[iter]] <- lapply(SSValid_win, '[[',iter)
    # compute spectral species distribution for each cluster & dissimilarity
    ssd_val <- lapply(SSdist,
                      FUN = get_beta_from_ssd,
                      beta_metrics = beta_metrics,
                      nb_clusters = nb_clusters,
                      pcelim = pcelim)
    
    mantelVal <- CorrVal <- Assess <- list()
    for (beta in beta_metrics){
      mat_diss_iter <- lapply(lapply(ssd_val, '[[','mat_diss'), '[[', beta)
      ssd <- lapply(ssd_val, '[[','ssd')
      mat_diss <- Reduce('+', mat_diss_iter)/nb_iter
      mat_diss_full <- matrix(data = NA,
                              nrow = nbPlots_total,
                              ncol = nbPlots_total)
      mat_diss_full[IDwindow,IDwindow] <- mat_diss
      colnames(mat_diss_full) <- rownames(mat_diss_full) <- Attributes$ID_biodivMapR
      # Corr_val <- cor.test(c(obs2optimize[['BC']]),c(MatBC_Full),
      #                      method = 'pearson')$estimate
      
      Assess[[beta]] <- mat_diss_full
      mantelVal[[beta]] <- vegan::mantel(xdis = as.dist(mat_diss_full),
                                 ydis = as.dist(obs2optimize[[beta]]),
                                 na.rm = TRUE, method = corrMethod)
      CorrVal[[beta]] <- mantelVal[[beta]]$statistic
      # CorrVal[['BC']]  <- cor.test(obs2optimize[['BC']],
      #                              c(Assess[['BC']]), method = 'pearson')$estimate
    }
  } else {
    Assess <- CorrVal  <- NA
  }
  
  # if (!is.null(obs2optimize$BC)){
  #   # compute BC matrix from spectral species
  #   SSValid_win <- SSValid %>% group_split(win_ID, .keep = FALSE)
  #   # spectral species distribution
  #   SSdist <- list()
  #   for (iter in names(SSValid_win[[1]]))
  #     SSdist[[iter]] <- lapply(SSValid_win, '[[',iter)
  #   # compute spectral species distribution for each cluster & BC dissimilarity
  #   SSD_BCval <- lapply(SSdist,
  #                       FUN = get_bc_diss_from_ssd,
  #                       nb_clusters = nb_clusters,
  #                       pcelim = pcelim)
  #   MatBC_iter <- lapply(SSD_BCval, '[[','MatBC')
  #   SSD <- lapply(SSD_BCval, '[[','SSD')
  #   MatBC <- Reduce('+', MatBC_iter)/nb_iter
  #   MatBC_Full <- matrix(data = NA,
  #                        nrow = nbPlots_total,
  #                        ncol = nbPlots_total)
  #   MatBC_Full[IDwindow,IDwindow] <- MatBC
  #   colnames(MatBC_Full) <- rownames(MatBC_Full) <- Attributes$ID_biodivMapR
  #   # Corr_val <- cor.test(c(obs2optimize[['BC']]),c(MatBC_Full),
  #   #                      method = 'pearson')$estimate
  # 
  #   Assess <- MatBC_Full
  #   mantelVal <- vegan::mantel(xdis = as.dist(MatBC_Full),
  #                              ydis = as.dist(obs2optimize[['BC']]),
  #                              na.rm = TRUE, method = corrMethod)
  #   CorrVal <- mantelVal$statistic
  #   # CorrVal[['BC']]  <- cor.test(obs2optimize[['BC']],
  #   #                              c(Assess[['BC']]), method = 'pearson')$estimate
  # } else {
  #   Assess <- CorrVal  <- NA
  # }
  return(list('Assess' = Assess, 'CorrVal' = CorrVal))
}
