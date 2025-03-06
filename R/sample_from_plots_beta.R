#' get samples for beta diversity mapping
#'
#' @param feature_dir character.
#' @param list_features character.
#' @param plots list.
#' @param nbPixValid numeric.
#' @param mask_dir character.
#' @param nbsamples_beta numeric.
#' @param window_size numeric.
#' @param nbCPU numeric.
#'
#' @return samples_beta_terra
#' @importFrom progressr progressor handlers with_progress
#' @export

sample_from_plots_beta <- function(feature_dir, list_features, plots,
                                   nbPixValid, mask_dir = NULL, window_size,
                                   nbsamples_beta = 2e3, nbCPU = 1){

  # define number of samples per tile
  totalPixels <- sum(unlist(nbPixValid))
  if (totalPixels<nbsamples_beta)
    nbsamples_beta <- totalPixels
  ratioStats <- nbsamples_beta/totalPixels
  plots2sel <- lapply(X = nbPixValid,
                      FUN = function(x, ratio){as.numeric(ceiling(x*ratio))},
                      ratio = ratioStats)

  # define features to sample
  if (is.null(mask_dir)){
    listfiles <- list.files(feature_dir, full.names = T)
    feat_list <- list_features
  } else {
    listfiles <- c(list.files(feature_dir, full.names = T),
                   list.files(mask_dir, full.names = T))
    feat_list <- c(list_features, 'mask')
  }

  # extract samples
  if (nbCPU==1){
    progressr::handlers("cli")
    suppressWarnings(with_progress({
      p <- progressr::progressor(steps = length(plots),
                                 message = 'get samples for beta diversity')
      samples_beta_terra <- mapply(FUN = get_plots_from_tiles,
                                   plotID = names(plots), plots2sel = plots2sel,
                                   MoreArgs = list(listfiles = listfiles,
                                                   feat_list = feat_list,
                                                   window_size = window_size,
                                                   p = p),
                                   SIMPLIFY = F)}))
  } else {
    message('get samples for beta diversity')
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)
    samples_beta_terra <- future.apply::future_mapply(FUN = get_plots_from_tiles,
                                                      plotID = names(plots), plots2sel = plots2sel,
                                                      MoreArgs = list(listfiles = listfiles,
                                                                      feat_list = feat_list,
                                                                      window_size = window_size),
                                                      future.seed = T,
                                                      future.chunk.size = NULL,
                                                      future.scheduling = structure(TRUE, ordering = "random"),
                                                      SIMPLIFY = F)
    parallel::stopCluster(cl)
    plan(sequential)
  }
  return(samples_beta_terra)
}
