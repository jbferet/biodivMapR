#' performs SFS to identify combination of input variables maximizing a
#' criterion
#'
#' @param input_raster SpatRaster or list of SpatRaster
#' @param obs_vect SpatVector or SpatVectorCollection
#' @param obs2optimize numeric .list of ground obs diversity metrics
#' corresponding to obs_vect.
#' Expected values: richness, shannon, simpson, hill, BC
#' @param obs_criterion character. richness, shannon, simpson or BC
#' @param corrMethod character. select between method available for cor.test
#' @param input_mask SpatRaster corresponding to mask
#' @param Hill_order numeric.
#' @param nb_clusters numeric.
#' @param min_sun numeric.
#' @param nb_pix numeric.
#' @param nb_iter numeric.
#' @param pcelim numeric.
#' @param verbose boolean.
#' @param nbWorkers numeric.
#' @param nbCPU numeric.
#'
#' @return list including performances (correlation) of SFS with additional
#' features and assessed diversity metrics corresponding to each step
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan multisession sequential
#' @importFrom foreach foreach %dopar%
#' @importFrom dplyr group_split
#' @importFrom vegan mantel
#' @importFrom progressr progressor handlers with_progress
#' @importFrom stats cor.test
#' @importFrom parallel makeCluster stopCluster
#'
#' @export

biodivMapR_sfs <- function(input_raster, obs_vect, obs2optimize,
                           obs_criterion = 'shannon', corrMethod = 'spearman',
                           input_mask = NULL, Hill_order = 1, nb_clusters = 50,
                           min_sun = 0.25, nb_pix = 1e5, nb_iter = 10,
                           pcelim = 0.02, verbose = TRUE, nbWorkers = 1,
                           nbCPU = 1){

  FullListIndices <- c('richness', 'shannon', 'simpson', 'hill', 'BC',
                       'FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')
  #### Which diversity metrics should be computed?
  alphamet <- c('richness', 'shannon', 'simpson', 'hill')
  betamet <- 'BC'
  fmet <- c('FRic', 'FEve', 'FDiv', 'FDis', 'FRaoq')
  # if computation of functional metrics required
  alpha_metrics <- alphamet[which(alphamet %in% names(obs2optimize))]
  if (length(alpha_metrics)==0)
    alpha_metrics <- NULL
  # computation of beta diversity required?
  betametrics <- betamet[which(betamet %in% names(obs2optimize))]
  if (length(betametrics)==0)
    getBeta <- FALSE
  if (length(betametrics)>0)
    getBeta <- TRUE
  functionalmetrics <- fmet[which(fmet %in% names(obs2optimize))]
  if (length(functionalmetrics)==0)
    fd_metrics <- NULL
  if (length(functionalmetrics)>0)
    fd_metrics <- functionalmetrics


  # 1- prepare for kmeans over the full spatial extent
  if (verbose)
    message('sampling pixels to compute spectral species')
  Pix_Per_Iter <- define_pixels_per_iter(input_rast = input_raster,
                                         input_mask = input_mask,
                                         nb_pix = nb_pix,
                                         nb_iter = nb_iter)
  extent_area <- get_raster_extent(input_raster[[1]])
  nb_samples <- Pix_Per_Iter*nb_iter
  rast_sample <- sample_from_raster(extent_area = extent_area,
                                    nb_samples = nb_samples,
                                    input_rast = input_raster,
                                    input_mask = input_mask)
  if (verbose)
    message('sampling done')

  # 2- extract information from field plots
  if (inherits(obs_vect, what = 'SpatVectorCollection')){
    # extract from list of SpatVectors in collection
    rastext <- lapply(X = obs_vect,
                      FUN = extract_vect_from_rast,
                      input_rast = input_raster,
                      input_mask = input_mask,
                      min_sun = min_sun, prog = FALSE)
    # update plot ID in collection
    nbPlots_total <- 0
    for (ind_vect in seq_len(length(obs_vect))){
      AttributeTable <- rastext[[ind_vect]]$AttributeTable
      rast_sample_vect <- rastext[[ind_vect]]$rast_sample_vect
      AttributeTable$ID_biodivMapR <- AttributeTable$ID_biodivMapR +
        nbPlots_total
      rast_sample_vect$ID <- rast_sample_vect$ID + nbPlots_total
      nbPlots_total <- max(AttributeTable$ID_biodivMapR)
      rastext[[ind_vect]]$AttributeTable <- AttributeTable
      rastext[[ind_vect]]$rast_sample_vect <- rast_sample_vect
    }
    rast_sample_vect <- lapply(rastext,'[[','rast_sample_vect')
    AttributeTable <- lapply(rastext,'[[','AttributeTable')
    rast_val <- do.call(rbind,rast_sample_vect)
    Attributes <- do.call(rbind,AttributeTable)
  } else if (inherits(obs_vect, what = 'SpatVector')){
    # extract from list of SpatVectors in collection
    rastext <- extract_vect_from_rast(SpatVector = obs_vect,
                                      input_rast = input_raster,
                                      input_mask = input_mask,
                                      min_sun = min_sun)
    # update plot ID in collection
    # rast_sample_vect <- lapply(rastext,'[[','rast_sample_vect')
    # AttributeTable <- lapply(rastext,'[[','AttributeTable')
    # rast_sample_vect <- rastext$rast_sample_vect
    # AttributeTable <- rastext$AttributeTable
    rast_val <- rastext$rast_sample_vect
    Attributes <- rastext$AttributeTable
    nbPlots_total <- length(Attributes$id)
  }
  if (verbose)
    message('plot extraction done')
  IDplot <- rast_val$ID
  rast_val$ID <- NULL

  # 3- perform SFS
  if (verbose)
    message('perform SFS')
  # component visual pre-selection
  # AllVars <- seq_len(ncol(rast_sample))
  AllVars <- names(rast_sample)
  nb_pcs_to_keep <- length(AllVars)
  # multi-thread feature selection (SFS)
  registerDoFuture()
  # plan(multisession, workers = nbWorkers)
  cl <- parallel::makeCluster(nbWorkers)
  with(plan('cluster', workers = cl), local = TRUE)

  Corr_criterion <- EvolCorr <-   CorrSFS <- AssessSFS <- list()
  SelectedVars <- EvolCorr$richness <- EvolCorr$shannon <- EvolCorr$simpson <-
    EvolCorr$hill <- EvolCorr$FRic <- EvolCorr$FEve <- EvolCorr$FDiv <-
    EvolCorr$BC <- c()

  progressr::handlers(global = TRUE)
  progressr::handlers('cli')
  progressr::with_progress({
    p <- progressr::progressor(steps = nb_pcs_to_keep)

  # pb <- progress::progress_bar$new(
  #   format = 'Perform feature selection [:bar] :percent in :elapsedfull',
  #   total = nb_pcs_to_keep, clear = FALSE, width= 100)

    for (nbvars2select in seq_len(nb_pcs_to_keep)){
      NumVar_list <- as.list(seq_len(length(AllVars)))
      numvar <- win_ID <- NULL
      subfeatures_SFS <- function() {
        foreach(numvar = NumVar_list) %dopar% {
          CorrVal <- Assess <- list()
          SelFeat_tmp <- c(SelectedVars,AllVars[[numvar]])
          # select features to compute kmeans
          Kmeans_info <- get_kmeans(rast_sample = rast_sample[SelFeat_tmp],
                                    nb_iter = nb_iter,
                                    nb_clusters = nb_clusters,
                                    nbCPU = 1, progressbar = FALSE)
          if (!is.null(fd_metrics)){
            # center reduce data
            inputdata_cr <- center_reduce(x = rast_val[SelFeat_tmp],
                                          m = Kmeans_info$MinVal,
                                          sig = Kmeans_info$Range)
            inputdata_cr$ID <- IDplot
            inputdata_cr <- inputdata_cr %>% split(.$ID)
            inputdata_cr <- lapply(inputdata_cr,
                                   function(x) data.frame(x[ , !(names(x) %in% 'ID')]))
            FunctDiv <- lapply(X = inputdata_cr,
                               FUN = get_functional_diversity,
                               fd_metrics = fd_metrics)
            FunctDiv <- data.frame('FRic' = unlist(lapply(FunctDiv, '[[',1)),
                                   'FEve' = unlist(lapply(FunctDiv, '[[',2)),
                                   'FDiv' = unlist(lapply(FunctDiv, '[[',3)))

            for (idx in fd_metrics){
              Assess[[idx]] <- FunctDiv[[idx]]
              CorrVal[[idx]]  <- cor.test(obs2optimize[[idx]],
                                          Assess[[idx]],
                                          method = corrMethod)$estimate
            }
          }
          # if (obs_criterion %in% c('richness', 'shannon', 'simpson', 'hill', 'BC')){
          SSValid <- get_spectralSpecies(inputdata = rast_val[SelFeat_tmp],
                                         Kmeans_info = Kmeans_info)
          SSValid$win_ID <- IDplot
          windows_per_plot <- split_chunk(SSchunk = SSValid, nbCPU = 1)
          windows_per_plot$win_ID <- list(SSValid$win_ID)
          alphabetaIdx_CPU <- lapply(X = windows_per_plot$SSwindow_perCPU,
                                     FUN = alphabeta_window_list,
                                     nb_clusters = nb_clusters,
                                     alpha_metrics = c('richness', 'shannon',
                                                       'simpson', 'hill'),
                                     Hill_order = Hill_order, pcelim = pcelim)
          alphabetaIdx <- unlist(alphabetaIdx_CPU,recursive = FALSE)
          rm(alphabetaIdx_CPU)
          gc()
          # 7- reshape alpha diversity metrics
          IDwindow <- unlist(windows_per_plot$IDwindow_perCPU)
          selcrit <- list('richness'= 1,
                          'shannon' = 3,
                          'simpson' = 5,
                          'hill' = 9)
          for (crit in names(selcrit)){
            if (!is.null(obs2optimize[[crit]])){
              Assess[[crit]] <- rep(x = NA,nbPlots_total)
              Assess[[crit]][IDwindow] <- unlist(lapply(alphabetaIdx,'[[',
                                                        selcrit[[crit]]))
              CorrVal[[crit]]  <- cor.test(obs2optimize[[crit]],
                                           Assess[[crit]],
                                           method = corrMethod)$estimate
            } else {
              CorrVal[[crit]] <- Assess[[crit]] <- NA
            }
          }
          if (!is.null(obs2optimize$BC)){
            # compute BC matrix from spectral species
            SSValid_win <- SSValid %>% group_split(win_ID, .keep = FALSE)
            # spectral species distribution
            SSdist <- list()
            for (iter in names(SSValid_win[[1]]))
              SSdist[[iter]] <- lapply(SSValid_win, '[[',iter)
            # compute spectral species distribution for each cluster & BC dissimilarity
            SSD_BCval <- lapply(SSdist,
                                FUN = get_bc_diss_from_ssd,
                                nb_clusters = nb_clusters,
                                pcelim = pcelim)
            MatBC_iter <- lapply(SSD_BCval, '[[','MatBC')
            SSD <- lapply(SSD_BCval, '[[','SSD')
            MatBC <- Reduce('+', MatBC_iter)/nb_iter
            MatBC_Full <- matrix(data = NA,
                                 nrow = nbPlots_total,
                                 ncol = nbPlots_total)
            MatBC_Full[IDwindow,IDwindow] <- MatBC
            colnames(MatBC_Full) <- rownames(MatBC_Full) <- Attributes$ID_biodivMapR
            # Corr_val <- cor.test(c(obs2optimize[['BC']]),c(MatBC_Full),
            #                      method = 'pearson')$estimate

            Assess[['BC']] <- MatBC_Full
            mantelVal <- vegan::mantel(xdis = as.dist(MatBC_Full),
                                       ydis = as.dist(obs2optimize[['BC']]),
                                       na.rm = TRUE, method = corrMethod)
            CorrVal[['BC']] <- mantelVal$statistic
            # CorrVal[['BC']]  <- cor.test(obs2optimize[['BC']],
            #                              c(Assess[['BC']]), method = 'pearson')$estimate
          } else {
            Assess[['BC']] <- CorrVal[['BC']]  <- NA
          }
          # }
          return(list('crit2Opt' = CorrVal,
                      'AssessedVal' = Assess))
        }
      }
      subSFS <- subfeatures_SFS()
      p(message = sprintf('Perform feature selection %g', nbvars2select))
      # pb$tick()
      CorrSFS[[nbvars2select]] <- list()
      for (ind in FullListIndices)
        CorrSFS[[nbvars2select]][[ind]] <- unlist(lapply(lapply(subSFS,
                                                                '[[','crit2Opt'),
                                                         '[[',ind))
      CorrSFS[[nbvars2select]] <- data.frame(CorrSFS[[nbvars2select]])
      rownames(CorrSFS[[nbvars2select]]) <- AllVars
      criterion <- CorrSFS[[nbvars2select]][[obs_criterion]]
      SelVar <- which(criterion == max(criterion,na.rm = T))
      AssessSFS[[nbvars2select]] <- list()
      for (ind in FullListIndices) {
        AssessSFS[[nbvars2select]][[ind]] <- unlist((lapply(lapply(subSFS,'[[',
                                                                   'AssessedVal'),
                                                            '[[',ind))[SelVar[1]])
      }
      # which criterion to maximize with SFS?
      WhichVar <- rownames(CorrSFS[[nbvars2select]])[SelVar[1]]
      # add selected component to selected vars
      SelectedVars <- c(SelectedVars,WhichVar)
      # delete selected component from AllVars
      AllVars <- AllVars[-which(AllVars==WhichVar)]
      EvolCorr$richness <- c(EvolCorr$richness,
                             CorrSFS[[nbvars2select]]$richness[SelVar[1]])
      EvolCorr$shannon <- c(EvolCorr$shannon,
                            CorrSFS[[nbvars2select]]$shannon[SelVar[1]])
      EvolCorr$simpson <- c(EvolCorr$simpson,
                            CorrSFS[[nbvars2select]]$simpson[SelVar[1]])
      EvolCorr$hill <- c(EvolCorr$hill,
                         CorrSFS[[nbvars2select]]$hill[SelVar[1]])
      EvolCorr$BC <- c(EvolCorr$BC,
                       CorrSFS[[nbvars2select]]$BC[SelVar[1]])

      EvolCorr$FRic <- c(EvolCorr$FRic, CorrSFS[[nbvars2select]]$FRic[SelVar[1]])
      EvolCorr$FEve <- c(EvolCorr$FEve, CorrSFS[[nbvars2select]]$FEve[SelVar[1]])
      EvolCorr$FDiv <- c(EvolCorr$FDiv, CorrSFS[[nbvars2select]]$FDiv[SelVar[1]])
    }
  })




  parallel::stopCluster(cl)
  plan(sequential)
  EvolCorr <- data.frame(EvolCorr)
  rownames(EvolCorr) <- SelectedVars
  AssessSFS_comp <- list()
  for (ind in FullListIndices) {
    if (!is.null(AssessSFS[[1]][[ind]])) {
      AssessSFS_comp[[ind]] <- data.frame(lapply(AssessSFS,'[[',ind))
      names(AssessSFS_comp[[ind]]) <- SelectedVars
      AssessSFS_comp[[ind]]$obs <- c(obs2optimize[[ind]])
    } else {
      AssessSFS_comp[[ind]] <- NA
    }
  }
  return(list('SFS_perf' = EvolCorr,
              'AssessDiv' = AssessSFS_comp))
}
