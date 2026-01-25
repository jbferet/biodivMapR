#' apply biodivMapR on a test set for different numbers of clusters
#'
#' @param input_raster SpatRaster or list of SpatRaster
#' @param obs_vect SpatVector or SpatVectorCollection
#' @param obs2optimize numeric .list of ground obs diversity metrics
#' corresponding to obs_vect.
#' Expected values: richness, shannon, simpson, BC
#' @param selected_bands numeric. bnds to select from input_raster
#' @param obs_criterion character. richness, shannon, simpson or BC
#' @param input_mask SpatRaster corresponding to mask
#' @param outputdir character. output directory
#' @param nb_clusters numeric.
#' @param min_sun numeric.
#' @param nb_iter numeric.
#' @param pcelim numeric.
#' @param verbose boolean.
#' @param nb_repetitions numeric.
#' @param nb_samples_alpha numeric.
#' @param Hill_order numeric.
#' @param algorithm character.
#' @param nbCPU numeric.
#'
#' @return list including performances (correlation) of SFS with additional
#' features and assessed diversity metrics corresponding to each step
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan multisession sequential
#' @importFrom foreach foreach %dopar%
#' @importFrom progressr progressor handlers with_progress
#' @importFrom stats cor.test
#' @importFrom parallel makeCluster stopCluster
#'
#' @export

biodivMapR_opt_clusters <- function(input_raster, obs_vect, obs2optimize,
                                    selected_bands, obs_criterion = 'shannon',
                                    input_mask = NULL, outputdir = './',
                                    nb_clusters = 50, min_sun = 0.25,
                                    nb_iter = 10, pcelim = 0.02, verbose = TRUE,
                                    nb_repetitions = 50, nb_samples_alpha = 1e5,
                                    Hill_order = 1, algorithm = 'Hartigan-Wong',
                                    nbCPU = 1){

  #### Which diversity metrics should be computed?
  alphamet <- c('richness', 'shannon', 'simpson', 'hill')
  betamet <- 'BC'
  fmet <- c('FRic', 'FEve', 'FDiv')
  # if computation of functional metrics required
  alpha_metrics <- alphamet[which(alphamet %in% obs_criterion)]
  if (length(alpha_metrics)==0)
    alpha_metrics <- NULL
  # computation of beta diversity required?
  betametrics <- betamet[which(betamet %in% obs_criterion)]
  if (length(betametrics)==0)
    getBeta <- FALSE
  if (length(betametrics)>0)
    getBeta <- TRUE
  # computation of functional metrics
  fd_metrics <- NULL

  # prepare sequence of clusters to test over multiple repetitions
  if (length(nb_clusters)==1)
    nbClust_list <- seq(2,nb_clusters)
  if (length(nb_clusters)>1)
    nbClust_list <- nb_clusters
  divIndex_est <- list()

  progressr::handlers(global = TRUE)
  progressr::handlers('cli')
  progressr::with_progress({
    p <- progressr::progressor(steps = nb_repetitions)
    # pb <- progress_bar$new(
    #   format = "Repeat clustering [:bar] :percent in :elapsedfull",
    #   total = nb_repetitions, clear = FALSE, width= 100)

    for (repet in seq_len(nb_repetitions)){
      # extract information from SpatVectorCollection or SpatVector
      if (inherits(obs_vect, what = 'SpatVectorCollection')){
        rastext <- extract_svc_from_rast(SpatVector = obs_vect,
                                         input_rast = input_raster,
                                         input_mask = input_mask,
                                         min_sun = min_sun, prog = FALSE)
        rast_val <- rastext$rast_sample_vect
        Attributes <- rastext$AttributeTable
        nbPlots_total <- nrow(Attributes)
      } else if (inherits(obs_vect, what = 'SpatVector')){
        rastext <- extract_vect_from_rast(SpatVector = obs_vect,
                                          input_rast = input_raster,
                                          input_mask = input_mask,
                                          min_sun = min_sun, prog = FALSE)
        rast_val <- rastext$rast_sample_vect
        Attributes <- rastext$AttributeTable
        nbPlots_total <- nrow(Attributes)
      }

      # perform kmeans for the full list of nbClust_list
      Kmeans_info <- explore_kmeans(input_rast = input_raster,
                                    input_mask = input_mask,
                                    selected_bands = selected_bands,
                                    nbClust_list = nbClust_list,
                                    nb_iter = nb_iter,
                                    nb_samples_alpha = nb_samples_alpha,
                                    algorithm = algorithm,
                                    nbCPU = nbCPU)
      gc()
      registerDoFuture()
      cl <- parallel::makeCluster(nbCPU)
      with(plan("cluster", workers = cl), local = TRUE)	
      kmit <- NULL
      get_diversity_from_plots_list <- function() {
        foreach(kmit = Kmeans_info) %dopar% {
          divplots <- get_diversity_from_plots(input_rast = input_raster,
                                               validation_vect = obs_vect,
                                               Hill_order = Hill_order,
                                               Kmeans_info = kmit,
                                               Beta_info = NULL,
                                               input_mask  = input_mask,
                                               alpha_metrics = alpha_metrics,
                                               fd_metrics = fd_metrics,
                                               getBeta = getBeta,
                                               rast_sample = rast_val,
                                               AttributeTable = Attributes,
                                               selected_bands = selected_bands,
                                               min_sun = min_sun,
                                               pcelim = pcelim,
                                               nbCPU = 1)
          return(divplots)
        }
      }
      divPlots_kmeans <- get_diversity_from_plots_list()
      plan(sequential)
      p(message = sprintf("Repeat clustering %g", repet))
      # pb$tick()
      gc()

      divIndex_est[[repet]] <- list()
      for (crit0 in obs_criterion){
        if (crit0 %in% c('richness', 'shannon', 'simpson', 'hill')){
          Sel1 <- lapply(X = divPlots_kmeans, FUN = '[[', 'specdiv')
          divIndex_est[[repet]][[crit0]] <- data.frame(lapply(X = Sel1,
                                                              FUN = '[[',
                                                              paste0(crit0, '_mean')))
          colnames(divIndex_est[[repet]][[crit0]]) <- unlist(nbClust_list)
        } else if (crit0 %in% 'BC'){
          divIndex_est[[repet]][[crit0]] <- lapply(X = divPlots_kmeans,
                                                   FUN = '[[', 'BC_dissimilarity')
          divIndex_est[[repet]][[crit0]] <- lapply(X = divIndex_est[[repet]][[crit0]],
                                                   FUN = 'as.dist')
          divIndex_est[[repet]][[crit0]] <- data.frame(lapply(X = divIndex_est[[repet]][[crit0]],
                                                              FUN = 'c'))
          colnames(divIndex_est[[repet]][[crit0]]) <- unlist(nbClust_list)
        }
        # save results to ease analysis if stopped at some point
        fileName <- file.path(outputdir, paste0(crit0,'_repet_',repet,'.csv'))
        readr::write_delim(x = round(divIndex_est[[repet]][[crit0]], digits = 5),
                           file = fileName, delim = '\t', progress = FALSE)
      }
    }
  })

  for (crit0 in obs_criterion){
    Est_Val_indiv <- lapply(X = divIndex_est, FUN = '[[', crit0)
    PearsonStats <- SpearmanStats <- data.frame('mean' = NULL, 'sd' = NULL)
    PearsonAll <- SpearmanAll <- list()
    for (nbc in unlist(nbClust_list)){
      # for each number of clusters, get all estimated values
      Est_Val_indiv2 <- lapply(X = Est_Val_indiv, FUN = '[[', as.character(nbc))
      # compute correlation with target
      corall <- lapply(X = Est_Val_indiv2, FUN = cor.test,
                       y = obs2optimize[[crit0]])
      Pearson <- unlist(lapply(X = corall, '[[', 'estimate'))
      PearsonAll[[as.character(nbc)]] <- Pearson
      PearsonStats <- rbind(PearsonStats,
                            data.frame('mean' = mean(Pearson),
                                       'sd' = sd(Pearson)))
      corall <- lapply(X = Est_Val_indiv2, FUN = cor.test,
                       y = obs2optimize[[crit0]],
                       method = 'spearman')
      Spearman <- unlist(lapply(X = corall, '[[', 'estimate'))
      SpearmanStats <- rbind(SpearmanStats, data.frame('mean' = mean(Spearman),
                                                       'sd' = sd(Spearman)))
      SpearmanAll[[as.character(nbc)]] <- Spearman
    }
    PearsonStats$nb_clusters <- unlist(nbClust_list)
    fileName <- file.path(outputdir, paste0(crit0, '_Pearson_mean.csv'))
    readr::write_delim(x = round(PearsonStats, digits = 5),
                       file = fileName, delim = '\t', progress = FALSE)

    SpearmanStats$nb_clusters <- unlist(nbClust_list)
    fileName <- file.path(outputdir, paste0(crit0, '_Spearman_Mean.csv'))
    readr::write_delim(x = round(SpearmanStats, digits = 5),
                       file = fileName, delim = '\t', progress = FALSE)

    SpearmanAll <- data.frame(do.call('rbind',SpearmanAll))
    fileName <- file.path(outputdir, paste0(crit0, '_Spearman_All.csv'))
    readr::write_delim(x = round(SpearmanAll, digits = 5),
                       file = fileName, delim = '\t', progress = FALSE)

    PearsonAll <- data.frame(do.call('rbind',PearsonAll))
    fileName <- file.path(outputdir, paste0(crit0, '_Pearson_All.csv'))
    readr::write_delim(x = round(PearsonAll, digits = 5),
                       file = fileName, delim = '\t', progress = FALSE)
  }
  return(list('Pearson' = PearsonStats,
              'Spearman' = SpearmanStats))
}
