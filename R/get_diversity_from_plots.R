#' computes diversity metrics from validation plots
#'
#' @param input_rast SpatRaster
#' @param validation_vect SpatVector
#' @param Hill_order numeric. Hill order
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics from training set
#' @param input_mask SpatRaster
#' @param Functional character.
#' @param SelectBands numeric. bands selected from input_rast
#' @param rast_sample dataframe.
#' @param AttributeTable dataframe.
#' @param alphametrics character.
#' @param MinSun numeric. minimum amount of sunlit pixels in the plots
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param getBeta boolean. set true if computation of beta required
#' @param verbose boolean. set true for messages
#'
#' @return SpatVector including diversity metrics and BC dissimilarity for the plots
#' @importFrom dplyr group_split
#' @importFrom stats as.dist
#' @export

get_diversity_from_plots <- function(input_rast, validation_vect,
                                     Hill_order = 1,
                                     Kmeans_info, Beta_info = NULL,
                                     input_mask  = NULL, Functional = NULL,
                                     SelectBands = NULL,
                                     rast_sample = NULL, AttributeTable = NULL,
                                     alphametrics = c('richness', 'shannon', 'simpson', 'hill'),
                                     MinSun = 0.25, pcelim = 0.02, nbCPU = 1, getBeta = T,
                                     verbose = F){
  if (verbose == T) message('Compute diversity from vector plot network')
  FunctDiv <- MatBC_Full <- win_ID <- NULL
  # get nbIter and nbclusters
  nbIter <- length(Kmeans_info$Centroids)
  nbclusters <- dim(Kmeans_info$Centroids[[1]])[1]
  # get dimPCO
  if (is.null(Beta_info)){
    dimPCO <- 3
  } else {
    dimPCO <- ncol(Beta_info$BetaPCO$points)
  }
  # read vector data
  if (inherits(validation_vect, what = 'SpatVectorCollection') & is.null(rast_sample)){
    SSValid <- Attributes <- list()
    FRic <- FEve <- FDiv <- FDis <- FRaoq <- list()
    nbPlots_init <- 0
    for (ind_vect in seq_len(length(validation_vect))){
      ssvect <- spectralspecies_per_polygon(SpatVector = validation_vect[[ind_vect]],
                                            input_rast = input_rast,
                                            Functional = Functional,
                                            SelectBands = SelectBands,
                                            input_mask = input_mask,
                                            Kmeans_info = Kmeans_info,
                                            rast_sample = rast_sample,
                                            AttributeTable = AttributeTable,
                                            MinSun = MinSun)
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
                                          Functional = Functional,
                                          Kmeans_info = Kmeans_info,
                                          SelectBands = SelectBands,
                                          rast_sample = rast_sample,
                                          AttributeTable = AttributeTable,
                                          MinSun = MinSun)
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

  windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
  windows_per_plot$win_ID <- list(SSValid$win_ID)

  alphabetaIdx_CPU <- lapply(X = windows_per_plot$SSwindow_perCPU,
                             FUN = alphabeta_window_list,
                             nbclusters = nbclusters,
                             alphametrics = alphametrics,
                             Hill_order = Hill_order,
                             Beta_info = Beta_info, pcelim = pcelim)

  alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = F)
  rm(alphabetaIdx_CPU)
  gc()
  # 7- reshape alpha diversity metrics
  IDwindow <- unlist(windows_per_plot$IDwindow_perCPU)
  res_shapeChunk <- list()
  for (i in seq_len(10)) {
    restmp <- unlist(lapply(alphabetaIdx,'[[',i))
    res_shapeChunk[[i]] <- rep(x = NA,nbPlots_init)
    res_shapeChunk[[i]][IDwindow] <- restmp
  }

  Attributes$richness_mean <- res_shapeChunk[[1]]
  Attributes$richness_sd <- res_shapeChunk[[2]]
  Attributes$shannon_mean <- res_shapeChunk[[3]]
  Attributes$shannon_sd <- res_shapeChunk[[4]]
  Attributes$simpson_mean <- res_shapeChunk[[5]]
  Attributes$simpson_sd <- res_shapeChunk[[6]]
  # Attributes$fisher_mean <- res_shapeChunk[[7]]
  # Attributes$fisher_sd <- res_shapeChunk[[8]]
  Attributes$hill_mean <- res_shapeChunk[[9]]
  Attributes$hill_sd <- res_shapeChunk[[10]]
  # 8- reshape beta diversity metrics
  if (!is.null(Beta_info)){
    PCoA_BC0 <- do.call(rbind,lapply(alphabetaIdx,'[[',11))
    PCoA_BC <- matrix(data = NA,nrow = nbPlots_init, ncol = dimPCO)
    PCoA_BC[IDwindow,] <- PCoA_BC0
    Attributes$BetaFull_PCoA_1 <- PCoA_BC[,1]
    Attributes$BetaFull_PCoA_2 <- PCoA_BC[,2]
    Attributes$BetaFull_PCoA_3 <- PCoA_BC[,3]
  }

  if (getBeta ==T){
    # compute BC matrix from spectral species
    SSValid_win <- SSValid %>% group_split(win_ID, .keep = F)
    # spectral species distribution
    SSdist <- list()
    for (iter in names(SSValid_win[[1]])) SSdist[[iter]] <- lapply(SSValid_win, '[[',iter)
    # compute spectral species distribution for each cluster & BC dissimilarity
    SSD_BCval <- lapply(SSdist,
                        FUN = get_BCdiss_from_SSD,
                        nbclusters = nbclusters,
                        pcelim = pcelim)

    MatBC_iter <- lapply(SSD_BCval, '[[','MatBC')
    SSD <- lapply(SSD_BCval, '[[','SSD')
    MatBC <- Reduce('+', MatBC_iter)/nbIter
    MatBC_Full <- matrix(data = NA, nrow = nbPlots_init, ncol = nbPlots_init)
    MatBC_Full[IDwindow,IDwindow] <- MatBC
    MatBCdist <- stats::as.dist(MatBC, diag = FALSE, upper = FALSE)
    colnames(MatBC_Full) <- rownames(MatBC_Full) <- Attributes$ID_biodivMapR
    BetaPCO <- pco(MatBCdist, k = dimPCO)
    PCoA_BC <- matrix(data = NA,nrow = nbPlots_init, ncol = dimPCO)
    PCoA_BC[IDwindow,] <- BetaPCO$points
    Attributes$BetaPlots_PCoA_1 <- PCoA_BC[,1]
    Attributes$BetaPlots_PCoA_2 <- PCoA_BC[,2]
    Attributes$BetaPlots_PCoA_3 <- PCoA_BC[,3]
  }
  if (!is.null(Functional)) {
    # FunctDiv$ID_biodivMapR <- Attributes$ID_biodivMapR
    # FunctDiv$id <- Attributes$id
    # FunctDiv$source <- Attributes$source
    Attributes$FRic <- Attributes$FEve <- Attributes$FDiv <- Attributes$FDis <- Attributes$FRaoq <- NA
    Attributes$FRic[selPlots] <- FunctDiv$FRic
    Attributes$FEve[selPlots] <- FunctDiv$FEve
    Attributes$FDiv[selPlots] <- FunctDiv$FDiv
    Attributes$FDis[selPlots] <- FunctDiv$FDis
    Attributes$FRaoq[selPlots] <- FunctDiv$FRaoq
  }
  if (verbose == T) message('diversity computed from vector plot network')
  return(list('specdiv' = Attributes,
              'BC_dissimilarity' = MatBC_Full))
}
