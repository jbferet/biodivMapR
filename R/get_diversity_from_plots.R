#' computes diversity metrics from validation plots
#'
#' @param input_rast SpatRaster
#' @param validation_vect SpatVector
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param input_mask SpatRaster
#' @param selected_bands numeric. bands selected from input_rast
#' @param rast_sample dataframe.
#' @param alpha_metrics character.
#' @param Hill_order numeric. Hill order
#' @param getBeta boolean. set true if computation of beta required
#' @param beta_metrics character.
#' @param fd_metrics character.
#' @param AttributeTable dataframe.
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param verbose boolean. set true for messages
#'
#' @return SpatVector including diversity metrics and BC dissimilarity for the plots
#' @importFrom dplyr group_split
#' @importFrom stats as.dist
#' @export

get_diversity_from_plots <- function(
    input_rast, validation_vect, Kmeans_info, Beta_info = NULL,
    input_mask  = NULL, selected_bands = NULL, rast_sample = NULL,
    alpha_metrics = c('richness', 'shannon', 'simpson', 'hill'),
    Hill_order = 1, getBeta = TRUE,
    beta_metrics = c('bray', 'brayturn', 'jaccard', 'jaccardturn', 'sorensen', 'simpson_diss'),
    fd_metrics = NULL, AttributeTable = NULL, min_sun = 0.25, pcelim = 0.02,
    nbCPU = 1, verbose = FALSE){

  if (verbose)
    message('Compute diversity from vector plot network')
  FunctDiv <- MatBC_Full <- win_ID <- NULL
  # get nb_iter and nb_clusters
  nb_iter <- length(Kmeans_info$Centroids)
  nb_clusters <- dim(Kmeans_info$Centroids[[1]])[1]
  # get dimPCO
  if (is.null(Beta_info)){
    dimPCO <- 3
  } else {
    dimPCO <- ncol(Beta_info$beta_pco[[1]]$points)
  }
  # read vector data
  if (inherits(validation_vect,
               what = 'SpatVectorCollection') & is.null(rast_sample)){
    SSValid <- Attributes <- list()
    FRic <- FEve <- FDiv <- FDis <- FRaoq <- list()
    nbPlots_init <- 0
    for (ind_vect in seq_len(length(validation_vect))){
      ssvect <- spectralspecies_per_polygon(SpatVector = validation_vect[[ind_vect]],
                                            input_rast = input_rast,
                                            fd_metrics = fd_metrics,
                                            selected_bands = selected_bands,
                                            input_mask = input_mask,
                                            Kmeans_info = Kmeans_info,
                                            rast_sample = rast_sample,
                                            AttributeTable = AttributeTable,
                                            min_sun = min_sun)
      if (!is.null(ssvect$SSValid)){
        SSValid[[ind_vect]] <- ssvect$SSValid
        Attributes[[ind_vect]] <- ssvect$AttributeTable
        SSValid[[ind_vect]]$win_ID <- SSValid[[ind_vect]]$win_ID + nbPlots_init
        Attributes[[ind_vect]]$ID_biodivMapR <- Attributes[[ind_vect]]$ID_biodivMapR + nbPlots_init
        nbPlots_init <- nbPlots_init + length(validation_vect[[ind_vect]])
        FRic[[ind_vect]] <- ssvect$FunctDiv$FRic
        FEve[[ind_vect]] <- ssvect$FunctDiv$FEve
        FDiv[[ind_vect]] <- ssvect$FunctDiv$FDiv
        FDis[[ind_vect]] <- ssvect$FunctDiv$FDis
        FRaoq[[ind_vect]] <- ssvect$FunctDiv$FRaoq
      }
    }
    SSValid <- do.call(rbind,SSValid)
    Attributes <- do.call(rbind,Attributes)
    FunctDiv <- data.frame('FRic' = unlist(FRic),
                           'FEve' = unlist(FEve),
                           'FDiv' = unlist(FDiv),
                           'FDis' = unlist(FDis),
                           'FRaoq' = unlist(FRaoq))
  } else if (inherits(validation_vect, what = 'SpatVector') | (!is.null(rast_sample))){
    ssvect <- spectralspecies_per_polygon(SpatVector = validation_vect,
                                          input_rast = input_rast,
                                          input_mask = input_mask,
                                          fd_metrics = fd_metrics,
                                          Kmeans_info = Kmeans_info,
                                          selected_bands = selected_bands,
                                          rast_sample = rast_sample,
                                          AttributeTable = AttributeTable,
                                          min_sun = min_sun)
    SSValid <- ssvect$SSValid
    if (inherits(validation_vect, what = 'SpatVector')) {
      nbPlots_init <- length(validation_vect)
      nbPlots <- nrow(ssvect$AttributeTable)
      selPlots <- ssvect$AttributeTable$ID_biodivMapR
    } else if (!is.null(rast_sample)) {
      nbPlots_init <- nbPlots <- length(unique(rast_sample$ID))
      selPlots <- seq_len(nbPlots_init)
    }
    Attributes0 <- ssvect$AttributeTable
    Attributes <- data.frame(matrix(NA, ncol = ncol(ssvect$AttributeTable),
                                    # nrow = nrow(ssvect$AttributeTable)))
                                    nrow = nbPlots_init))
    names(Attributes) <- names(Attributes0)
    Attributes[selPlots,] <- Attributes0
    FunctDiv <- data.frame('FRic' = ssvect$FunctDiv$FRic,
                           'FEve' = ssvect$FunctDiv$FEve,
                           'FDiv' = ssvect$FunctDiv$FDiv,
                           'FDis' = ssvect$FunctDiv$FDis,
                           'FRaoq' = ssvect$FunctDiv$FRaoq)
  }
  # Attributes <- as.list(Attributes)

  windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
  windows_per_plot$win_ID <- list(SSValid$win_ID)

  alphabetaIdx <- lapply(X = windows_per_plot$SSwindow_perCPU,
                               # FUN = alphabeta_window_list,
                               FUN = alphabeta_chunk,
                             nb_clusters = nb_clusters,
                             alpha_metrics = alpha_metrics,
                             beta_metrics = beta_metrics,
                             Hill_order = Hill_order,
                             Beta_info = Beta_info, pcelim = pcelim, 
                             diss_output = TRUE)

  # alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = FALSE)
  # rm(alphabetaIdx_CPU)
  gc()
  # 7- reshape alpha diversity metrics
  IDwindow <- unlist(windows_per_plot$IDwindow_perCPU)
  
  names_alpha <- c('richness_mean', 'richness_sd', 'shannon_mean', 'shannon_sd', 
                   'simpson_mean', 'simpson_sd', 'hill_mean', 'hill_sd')
  res_shapeChunk <- list()
  for (aind in names_alpha) {
    restmp <- unlist(lapply(alphabetaIdx,'[[',aind))
    res_shapeChunk[[aind]] <- rep(x = NA,nbPlots_init)
    res_shapeChunk[[aind]][IDwindow] <- restmp
    Attributes[[aind]] <- res_shapeChunk[[aind]]
  }
  # Attributes$richness_mean <- res_shapeChunk[[1]]
  # Attributes$richness_sd <- res_shapeChunk[[2]]
  # Attributes$shannon_mean <- res_shapeChunk[[3]]
  # Attributes$shannon_sd <- res_shapeChunk[[4]]
  # Attributes$simpson_mean <- res_shapeChunk[[5]]
  # Attributes$simpson_sd <- res_shapeChunk[[6]]
  # Attributes$hill_mean <- res_shapeChunk[[9]]
  # Attributes$hill_sd <- res_shapeChunk[[10]]
  # 8- reshape beta diversity metrics
  names_beta <- c('pcoa_bray', 'pcoa_brayturn', 'pcoa_simpson_diss', 'pcoa_jaccard', 
                  'pcoa_jaccardturn', 'pcoa_sorensen')   
  if (!is.null(Beta_info)){
    for (bind in names_beta) {
      PCoA_BC0 <- do.call(rbind,lapply(alphabetaIdx,'[[',bind))
      if (!is.null(PCoA_BC0)){
        PCoA_BC <- matrix(data = NA,nrow = nbPlots_init, ncol = dimPCO)
        PCoA_BC[IDwindow,] <- PCoA_BC0
        Attributes[[paste0(bind, '_full_pcoa_1')]] <- PCoA_BC[,1]
        Attributes[[paste0(bind, '_full_pcoa_2')]] <- PCoA_BC[,2]
        Attributes[[paste0(bind, '_full_pcoa_3')]] <- PCoA_BC[,3]
      }
    }
  }

  if (getBeta){
    # compute BC matrix from spectral species
    SSValid_win <- SSValid %>% group_split(win_ID, .keep = FALSE)
    # spectral species distribution
    SSdist <- list()
    for (iter in names(SSValid_win[[1]]))
      SSdist[[iter]] <- lapply(SSValid_win, '[[',iter)
    # compute spectral species distribution for each cluster & BC dissimilarity
    SSD_BCval <- lapply(SSdist,
                        FUN = get_beta_from_ssd,
                        beta_metrics = beta_metrics, 
                        nb_clusters = nb_clusters,
                        pcelim = pcelim)

    dissimilarity_matrices <- list()
    for (beta in beta_metrics){
      MatBC_iter <- lapply(lapply(SSD_BCval, '[[','mat_diss'), '[[',beta)
      MatBC <- Reduce('+', MatBC_iter)/nb_iter
      dissimilarity_matrices[[beta]] <- matrix(data = NA, nrow = nbPlots_init, ncol = nbPlots_init)
      dissimilarity_matrices[[beta]][IDwindow,IDwindow] <- MatBC
      MatBCdist <- stats::as.dist(MatBC, diag = FALSE, upper = FALSE)
      colnames(dissimilarity_matrices[[beta]]) <- 
        rownames(dissimilarity_matrices[[beta]]) <- Attributes$ID_biodivMapR
      BetaPCO <- pco(MatBCdist, k = dimPCO)
      PCoA_BC <- matrix(data = NA,nrow = nbPlots_init, ncol = dimPCO)
      PCoA_BC[IDwindow,] <- BetaPCO$points
      Attributes[[paste0(beta, '_plot_pcoa_1')]] <- PCoA_BC[,1]
      Attributes[[paste0(beta, '_plot_pcoa_2')]] <- PCoA_BC[,2]
      Attributes[[paste0(beta, '_plot_pcoa_3')]] <- PCoA_BC[,3]
      
    }
  }
  if (!is.null(fd_metrics)) {
    # FunctDiv$ID_biodivMapR <- Attributes$ID_biodivMapR
    # FunctDiv$id <- Attributes$id
    # FunctDiv$source <- Attributes$source
    Attributes$FRic <- Attributes$FEve <- Attributes$FDiv <- NA
    Attributes$FDis <- Attributes$FRaoq <- NA
    Attributes$FRic[selPlots] <- FunctDiv$FRic
    Attributes$FEve[selPlots] <- FunctDiv$FEve
    Attributes$FDiv[selPlots] <- FunctDiv$FDiv
    Attributes$FDis[selPlots] <- FunctDiv$FDis
    Attributes$FRaoq[selPlots] <- FunctDiv$FRaoq
  }
  if (verbose)
    message('diversity computed from vector plot network')
  return(list('specdiv' = Attributes,
              'dissimilarity_matrices' = dissimilarity_matrices))
}
