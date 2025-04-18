#' apply biodivMapR on a test set for different numbers of clusters
#'
#' @param input_raster SpatRaster or list of SpatRaster
#' @param obs_vect SpatVector or SpatVectorCollection
#' @param obs2optimize numeric .list of ground obs diversity metrics
#' corresponding to obs_vect.
#' Expected values: richness, shannon, simpson, BC
#' @param SelectBands numeric. bnds to select from input_raster
#' @param obs_criterion character. richness, shannon, simpson or BC
#' @param input_mask SpatRaster corresponding to mask
#' @param outputdir character. output directory
#' @param nbclusters numeric.
#' @param MinSun numeric.
#' @param nbIter numeric.
#' @param pcelim numeric.
#' @param verbose boolean.
#' @param nbRep numeric.
#' @param maxPixel_kmeans numeric.
#' @param Hill_order numeric.
#' @param algorithm character.
#' @param nbCPU numeric.
#'
#' @return list including performances (correlation) of SFS with additional
#' features and assessed diversity metrics corresponding to each step
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan multisession sequential
#' @importFrom foreach foreach %dopar%
#' @importFrom vegan mantel
#' @importFrom progress progress_bar
#' @importFrom stats cor.test
#' @importFrom parallel makeCluster stopCluster
#'
#' @export

biodivMapR_OptClusters <- function(input_raster, obs_vect, obs2optimize, SelectBands,
                                   obs_criterion = 'shannon', input_mask = NULL,
                                   outputdir = './', nbclusters = 50,
                                   MinSun = 0.25, nbIter = 20,
                                   pcelim = 0.02, verbose = T, nbRep = 50,
                                   maxPixel_kmeans = 1e5,
                                   Hill_order = 1, algorithm = 'Hartigan-Wong',
                                   nbCPU = 1){

  #### Which diversity metrics should be computed?
  alphamet <- c('richness', 'shannon', 'simpson', 'hill')
  betamet <- 'BC'
  fmet <- c('FRic', 'FEve', 'FDiv')
  # if computation of functional metrics required
  alphametrics <- alphamet[which(alphamet %in% obs_criterion)]
  if (length(alphametrics)==0) alphametrics <- NULL
  # computation of beta diversity required?
  betametrics <- betamet[which(betamet %in% obs_criterion)]
  if (length(betametrics)==0) getBeta <- F
  if (length(betametrics)>0) getBeta <- T
  # computation of functional metrics required? No as no influence of nb clusters
  Functional <- NULL

  # prepare sequence of clusters to test over multiple repetitions
  if (length(nbclusters)==1) nbClust_list <- seq(2,nbclusters)
  if (length(nbclusters)>1) nbClust_list <- nbclusters
  divIndex_est <- list()
  pb <- progress_bar$new(
    format = "Repeat clustering [:bar] :percent in :elapsedfull",
    total = nbRep, clear = FALSE, width= 100)

  for (repet in seq_len(nbRep)){
    # extract information from SpatVectorCollection or SpatVector
    if (inherits(obs_vect, what = 'SpatVectorCollection')){
      rastext <- extract_svc_from_rast(SpatVector = obs_vect,
                                       input_rast = input_raster,
                                       input_mask = input_mask,
                                       MinSun = MinSun, prog = F)
      rast_val <- rastext$rast_sample_vect
      Attributes <- rastext$AttributeTable
      nbPlots_total <- nrow(Attributes)
    } else if (inherits(obs_vect, what = 'SpatVector')){
      rastext <- extract_vect_from_rast(SpatVector = obs_vect,
                                        input_rast = input_raster,
                                        input_mask = input_mask,
                                        MinSun = MinSun, prog = F)
      rast_val <- rastext$rast_sample_vect
      Attributes <- rastext$AttributeTable
      nbPlots_total <- nrow(Attributes)
    }

    # perform kmeans for the full list of nbClust_list
    Kmeans_info <- explore_kmeans(input_rast = input_raster,
                                  input_mask = input_mask,
                                  SelectBands = SelectBands,
                                  nbClust_list = nbClust_list,
                                  nbIter = nbIter,
                                  maxPixel_kmeans = maxPixel_kmeans,
                                  algorithm = algorithm,
                                  nbCPU = nbCPU)
    gc()
    registerDoFuture()
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)
	kmit <- NULL
    get_diversity_from_plots_list <- function() {
      foreach(kmit = Kmeans_info) %dopar% {
        divplots <- get_diversity_from_plots(input_rast = input_raster,
                                             validation_vect = obs_vect,
                                             Hill_order = Hill_order,
                                             Kmeans_info = kmit,
                                             Beta_info = NULL,
                                             input_mask  = input_mask,
                                             alphametrics = alphametrics,
                                             Functional = Functional,
                                             getBeta = getBeta,
                                             rast_sample = rast_val,
                                             AttributeTable = Attributes,
                                             SelectBands = SelectBands,
                                             MinSun = MinSun,
                                             pcelim = pcelim,
                                             nbCPU = 1)
        return(divplots)
      }
    }
    divPlots_kmeans <- get_diversity_from_plots_list()
    plan(sequential)
    pb$tick()
    gc()

    divIndex_est[[repet]] <- list()
    for (crit0 in obs_criterion){
      if (crit0 %in% c('richness', 'shannon', 'simpson', 'hill')){
        Sel1 <- lapply(X = divPlots_kmeans, FUN = '[[', 'validation_AlphaBeta')
        divIndex_est[[repet]][[crit0]] <- data.frame(lapply(X = Sel1, FUN = '[[', paste0(crit0, '_mean')))
        colnames(divIndex_est[[repet]][[crit0]]) <- unlist(nbClust_list)
      } else if (crit0 %in% 'BC'){
        divIndex_est[[repet]][[crit0]] <- lapply(X = divPlots_kmeans, FUN = '[[', 'BC_dissimilarity')
        divIndex_est[[repet]][[crit0]] <- lapply(X = divIndex_est[[repet]][[crit0]], FUN = 'as.dist')
        divIndex_est[[repet]][[crit0]] <- data.frame(lapply(X = divIndex_est[[repet]][[crit0]], FUN = 'c'))
        colnames(divIndex_est[[repet]][[crit0]]) <- unlist(nbClust_list)
      }
      # save results to ease analysis if stopped at some point
      fileName <- file.path(outputdir, paste0(crit0,'_repet_',repet,'.csv'))
      readr::write_delim(x = round(divIndex_est[[repet]][[crit0]], digits = 5),
                         file = fileName, delim = '\t', progress = F)
    }
  }
  for (crit0 in obs_criterion){
    Est_Val_indiv <- lapply(X = divIndex_est, FUN = '[[', crit0)
    PearsonStats <- SpearmanStats <- data.frame('mean' = NULL, 'sd' = NULL)
    PearsonAll <- SpearmanAll <- list()
    for (nbc in unlist(nbClust_list)){
      # for each number of clusters, get all estimated values
      Est_Val_indiv2 <- lapply(X = Est_Val_indiv, FUN = '[[', as.character(nbc))
      # compute correlation with target
      corall <- lapply(X = Est_Val_indiv2, FUN = cor.test, y = obs2optimize[[crit0]])
      Pearson <- unlist(lapply(X = corall, '[[', 'estimate'))
      PearsonAll[[as.character(nbc)]] <- Pearson
      PearsonStats <- rbind(PearsonStats, data.frame('mean' = mean(Pearson), 'sd' = sd(Pearson)))
      corall <- lapply(X = Est_Val_indiv2, FUN = cor.test, y = obs2optimize[[crit0]], method = 'spearman')
      Spearman <- unlist(lapply(X = corall, '[[', 'estimate'))
      SpearmanStats <- rbind(SpearmanStats, data.frame('mean' = mean(Spearman), 'sd' = sd(Spearman)))
      SpearmanAll[[as.character(nbc)]] <- Spearman
    }
    PearsonStats$nbClusters <- unlist(nbClust_list)
    fileName <- file.path(outputdir, paste0(crit0, '_Pearson_mean.csv'))
    readr::write_delim(x = round(PearsonStats, digits = 5),
                       file = fileName, delim = '\t', progress = F)

    SpearmanStats$nbClusters <- unlist(nbClust_list)
    fileName <- file.path(outputdir, paste0(crit0, '_Spearman_Mean.csv'))
    readr::write_delim(x = round(SpearmanStats, digits = 5),
                       file = fileName, delim = '\t', progress = F)

    SpearmanAll <- data.frame(do.call('rbind',SpearmanAll))
    fileName <- file.path(outputdir, paste0(crit0, '_Spearman_All.csv'))
    readr::write_delim(x = round(SpearmanAll, digits = 5),
                       file = fileName, delim = '\t', progress = F)

    PearsonAll <- data.frame(do.call('rbind',PearsonAll))
    fileName <- file.path(outputdir, paste0(crit0, '_Pearson_All.csv'))
    readr::write_delim(x = round(PearsonAll, digits = 5),
                       file = fileName, delim = '\t', progress = F)
  }
  return(list('Pearson' = PearsonStats,
              'Spearman' = SpearmanStats))
}
