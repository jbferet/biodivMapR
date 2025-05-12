#' get samples for alpha diversity mapping
#'
#' @param feature_dir character. directory where features to be used by biodivMapR are stored.
#' @param list_features list.
#' @param plots list.
#' @param nbCPU numeric.
#' @param nb_pix_valid list.
#' @param mask_dir character.
#' @param nb_samples_alpha numeric.
#' @param method character. method for terra::spatSample ('random' or 'regular')
#'
#' @return samples_alpha_terra
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_mapply
#' @export

sample_from_plots_alpha <- function(feature_dir, list_features, plots, nbCPU = 1,
                                    nb_pix_valid, mask_dir = NULL,
                                    nb_samples_alpha = 1e5, method = 'regular'){

  # define number of samples per tile
  totalPixels <- sum(unlist(nb_pix_valid))
  if (totalPixels<nb_samples_alpha)
    nb_samples_alpha <- totalPixels
  ratioStats <- nb_samples_alpha/totalPixels
  pix2sel <- lapply(X = nb_pix_valid,
                    FUN = function(x, ratio){as.numeric(round(x*ratio))},
                    ratio = ratioStats)

  # define features to sample
  if (is.null(mask_dir)){
    listfiles <- list.files(feature_dir, full.names = TRUE)
    listfiles <- unique(gsub(pattern = '.aux.xml', x = listfiles,
                             replacement = ''))
    feat_list <- list_features
  } else {
    listfiles1 <- list.files(feature_dir, full.names = TRUE)
    listfiles1 <- unique(gsub(pattern = '.aux.xml', x = listfiles1,
                              replacement = ''))
    listfiles2 <- list.files(mask_dir, full.names = TRUE)
    listfiles2 <- unique(gsub(pattern = '.aux.xml', x = listfiles2,
                              replacement = ''))
    listfiles <- c(listfiles1, listfiles2)
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
                                       method = method,
                                       as.df = TRUE, xy = FALSE,
                                       p = p),
                       SIMPLIFY = FALSE)}))
  } else {
    nbCPU2 <- min(c(4, nbCPU))
    message('get samples for alpha diversity')
    cl <- parallel::makeCluster(nbCPU2)
    plan("cluster", workers = cl)
    selpix <- future.apply::future_mapply(FUN = get_samples_from_tiles,
                                          plotID = names(plots),
                                          pix2sel = pix2sel,
                                          MoreArgs = list(listfiles = listfiles,
                                                          feat_list = feat_list,
                                                          method = method,
                                                          as.df = TRUE, xy = FALSE),
                                          future.seed = TRUE,
                                          future.chunk.size = NULL,
                                          future.scheduling = structure(TRUE, ordering = "random"),
                                          SIMPLIFY = FALSE)
    parallel::stopCluster(cl)
    plan(sequential)
  }
  # get alpha diversity metrics
  samples_alpha_terra <- do.call(what = 'rbind', selpix)
  samples_alpha_terra$ID <- seq_len(length(samples_alpha_terra[[1]]))
  return(samples_alpha_terra)
}
