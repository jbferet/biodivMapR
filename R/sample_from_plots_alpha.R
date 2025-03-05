#' get samples for alpha diversity mapping
#'
#' @param feature_dir character. directory where features to be used by biodivMapR are stored.
#' @param list_features list.
#' @param plots list.
#' @param nbCPU numeric.
#' @param nbPixValid list.
#' @param mask_dir character.
#' @param nbsamples_alpha numeric.
#'
#' @return samples_alpha_terra
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_mapply
#' @export

sample_from_plots_alpha <- function(feature_dir, list_features, plots, nbCPU = 1,
                                    nbPixValid, mask_dir = NULL,
                                    nbsamples_alpha = 1e5){

  # define number of samples per tile
  totalPixels <- sum(unlist(nbPixValid))
  if (totalPixels<nbsamples_alpha)
    nbsamples_alpha <- totalPixels
  ratioStats <- nbsamples_alpha/totalPixels
  pix2sel <- lapply(X = nbPixValid,
                    FUN = function(x, ratio){as.numeric(round(x*ratio))},
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
    handlers("cli")
    suppressWarnings(with_progress({
      p <- progressr::progressor(steps = length(plots),
                                 message = 'get samples for alpha diversity')
      selpix <- mapply(FUN = get_samples_from_tiles,
                       plotID = names(plots), pix2sel = pix2sel,
                       MoreArgs = list(listfiles = listfiles,
                                       feat_list = feat_list,
                                       as.df = T, xy = F,
                                       p = p),
                       SIMPLIFY = F)}))
  } else {
    message('get samples for alpha diversity')
    cl <- parallel::makeCluster(nbCPU)
    plan("cluster", workers = cl)
    selpix <- future.apply::future_mapply(FUN = get_samples_from_tiles,
                                          plotID = names(plots), pix2sel = pix2sel,
                                          MoreArgs = list(listfiles = listfiles,
                                                          feat_list = feat_list,
                                                          as.df = T, xy = F),
                                          future.seed = T, SIMPLIFY = F)
    parallel::stopCluster(cl)
    plan(sequential)
  }
  # get alpha diversity metrics
  samples_alpha_terra <- do.call(what = 'rbind', selpix)
  samples_alpha_terra$ID <- seq_len(length(samples_alpha_terra[[1]]))
  return(samples_alpha_terra)
}
