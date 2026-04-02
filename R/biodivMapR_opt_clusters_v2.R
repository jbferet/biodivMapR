#' apply biodivMapR on a test set for different numbers of clusters
#'
#' @param input_raster_path SpatRaster or list of SpatRaster
#' @param input_vector_path SpatVector or SpatVectorCollection
#' @param obs2optimize numeric .list of ground obs diversity metrics
#' corresponding to input_vector.
#' Expected values: richness, shannon, simpson, BC
#' @param selected_bands numeric. bnds to select from input_raster
#' @param obs_criterion character. richness, shannon, simpson or BC
#' @param input_mask_path SpatRaster corresponding to mask
#' @param outputdir character. output directory
#' @param nb_clusters numeric.
#' @param min_sun numeric.
#' @param nb_iter numeric.
#' @param pcelim numeric.
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
#' @importFrom utils write.table
#'
#' @export

biodivMapR_opt_clusters_v2 <- function(input_raster_path, input_vector_path, obs2optimize,
                                       selected_bands, obs_criterion = 'shannon',
                                       input_mask_path = NULL, outputdir = './',
                                       nb_clusters = 50, min_sun = 0.25,
                                       nb_iter = 10, pcelim = 0.02,
                                       nb_repetitions = 50, nb_samples_alpha = 1e5,
                                       Hill_order = 1, algorithm = 'Hartigan-Wong',
                                       nbCPU = 1){

  if (nbCPU > nb_repetitions)
    nbCPU <-  nb_repetitions

  #### Which diversity metrics should be computed?
  albet <- which_alpha_beta(obs_criterion)
  alpha_metrics <- albet$alpha_metrics
  # getBeta <- albet$getBeta

  # prepare sequence of clusters to test over multiple repetitions
  nbClust_list <- get_cluster_list(nb_clusters)
  list_repetitions <- as.list(num2char(seq_len(nb_repetitions),
                                       nbdigits = nchar(nb_repetitions)))

  # estimate diversity metrics over range of clusters number for each repetition
  if (nbCPU>1) {
    cl <- parallel::makeCluster(nbCPU)
    with(plan("cluster", workers = cl), local = TRUE)
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nb_repetitions)
      diversity_est <- future.apply::future_lapply(X = list_repetitions,
                                                   FUN = explore_cluster_range,
                                                   input_raster_path = input_raster_path,
                                                   input_vector_path = input_vector_path,
                                                   input_mask_path, selected_bands,
                                                   obs_criterion = obs_criterion,
                                                   alpha_metrics = alpha_metrics,
                                                   outputdir = outputdir,
                                                   nbClust_list = nbClust_list,
                                                   min_sun = min_sun,
                                                   nb_iter = nb_iter,
                                                   pcelim = pcelim,
                                                   nb_samples_alpha = nb_samples_alpha,
                                                   Hill_order = Hill_order,
                                                   algorithm = algorithm, p = p,
                                                   future.seed = TRUE)
    })
    parallel::stopCluster(cl)
    plan(sequential)
  } else {
    handlers("cli")
    with_progress({
      p <- progressr::progressor(steps = nb_repetitions)
      diversity_est <- lapply(X = list_repetitions,
                              FUN = explore_cluster_range,
                              input_raster_path = input_raster_path,
                              input_vector_path = input_vector_path,
                              input_mask_path, selected_bands,
                              obs_criterion = obs_criterion,
                              outputdir = outputdir,
                              nbClust_list = nbClust_list,
                              min_sun = min_sun,
                              nb_iter = nb_iter,
                              pcelim = pcelim,
                              nb_samples_alpha = nb_samples_alpha,
                              Hill_order = Hill_order,
                              algorithm = algorithm, p = p)
    })
  }

  pearson_stats <- spearman_stats <- list()
  for (crit0 in obs_criterion){
    Est_Val_indiv <- lapply(X = diversity_est, FUN = '[[', crit0)
    pearson_stats[[crit0]] <- spearman_stats[[crit0]] <- data.frame('mean' = NULL, 'sd' = NULL)
    pearson_all <- spearman_all <- list()
    for (nbc in nbClust_list){
      # for each number of clusters, get all estimated values
      Est_Val_indiv2 <- lapply(X = Est_Val_indiv, FUN = '[[', as.character(nbc))
      # compute correlation with target
      corall <- lapply(X = Est_Val_indiv2, FUN = cor.test,
                       y = obs2optimize[[crit0]])
      Pearson <- unlist(lapply(X = corall, '[[', 'estimate'))
      pearson_all[[as.character(nbc)]] <- Pearson
      pearson_stats[[crit0]] <- rbind(pearson_stats[[crit0]],
                                      data.frame('mean' = mean(Pearson),
                                                 'sd' = sd(Pearson)))
      corall <- lapply(X = Est_Val_indiv2, FUN = cor.test,
                       y = obs2optimize[[crit0]],
                       method = 'spearman')
      Spearman <- unlist(lapply(X = corall, '[[', 'estimate'))
      spearman_stats[[crit0]] <- rbind(spearman_stats[[crit0]],
                                       data.frame('mean' = mean(Spearman),
                                                  'sd' = sd(Spearman)))
      spearman_all[[as.character(nbc)]] <- Spearman
    }
    pearson_stats[[crit0]]$nb_clusters <- spearman_stats[[crit0]]$nb_clusters <- unlist(nbClust_list)
    pearson_stats[[crit0]]$metric <- spearman_stats[[crit0]]$metric <- crit0
    fileName <- file.path(outputdir, paste0(crit0, '_pearson_mean.csv'))
    write.table(file = fileName, x = format(pearson_stats[[crit0]], digits = 5),
                quote = F, row.names = F, col.names = T, sep = '\t')

    fileName <- file.path(outputdir, paste0(crit0, '_spearman_mean.csv'))
    write.table(file = fileName, x = format(spearman_stats[[crit0]], digits = 5),
                quote = F, row.names = F, col.names = T, sep = '\t')

    pearson_all <- data.frame(do.call('rbind',pearson_all))
    fileName <- file.path(outputdir, paste0(crit0, '_pearson_all.csv'))
    write.table(file = fileName, x = format(pearson_all, digits = 5),
                quote = F, row.names = F, col.names = T, sep = '\t')

    spearman_all <- data.frame(do.call('rbind',spearman_all))
    fileName <- file.path(outputdir, paste0(crit0, '_spearman_all.csv'))
    write.table(file = fileName, x = format(spearman_all, digits = 5),
                quote = F, row.names = F, col.names = T, sep = '\t')
  }
  list_out <- list('pearson' = pearson_stats,
                   'spearman' = spearman_stats)
  rm(list = ls(all=TRUE)[-which(ls()=='list_out')])
  return(list('pearson' = pearson_stats,
              'spearman' = spearman_stats))
}
