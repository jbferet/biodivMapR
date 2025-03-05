#' computes mask for a series of rasters based on IQR
#'
#' @param feature_dir character.
#' @param feature_list character.
#' @param mask_dir character.
#' @param plots list.
#' @param weightIRQ numeric.
#' @param filetype character.
#' @param nbCPU numeric.
#' @param nbPixstats numeric.
#'
#' @return mask_path
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom parallel makeCluster stopCluster
#' @export
#'
compute_mask_IQR_tiles <- function(feature_dir, feature_list, mask_dir, plots,
                                   weightIRQ = 4, filetype = 'GTiff', nbCPU = 1,
                                   nbPixstats = 5e6){
  
  # first check if mask files exist
  mask_path <- file.path(mask_dir, paste0('mask_', names(plots), '_IQR.tiff'))
  which_process <- which(!file.exists(mask_path))
  if (length(which_process)>0){
    plots <- plots[which_process]
    # list files
    listfiles <- list.files(feature_dir, full.names = T)
    ##############################################################################
    # get number of pixels per tile
    if (nbCPU==1){
      handlers("cli")
      suppressWarnings(with_progress({
        p <- progressr::progressor(steps = length(plots),
                                   message = 'get valid pixels from tiles')
        nbPixValid <- lapply(X = names(plots), FUN = get_valid_pixels,
                             listfiles = listfiles, p = p)}))
    } else {
      message('get valid pixels from tiles')
      cl <- parallel::makeCluster(nbCPU)
      plan("cluster", workers = cl)
      nbPixValid <- future.apply::future_lapply(X = names(plots),
                                                FUN = get_valid_pixels,
                                                listfiles = listfiles,
                                                future.seed = TRUE)
      parallel::stopCluster(cl)
      plan(sequential)
    }
    # get total number of pixels
    totalPixels <- sum(unlist(nbPixValid))
    
    ##############################################################################
    # get statistics on nbPixstats
    if (totalPixels<nbPixstats)
      nbPixstats <- totalPixels
    ratioStats <- nbPixstats/totalPixels
    pix2sel <- lapply(X = nbPixValid,
                      FUN = function(x, ratio){as.numeric(round(x*ratio))},
                      ratio = ratioStats)
    if (nbCPU==1){
      handlers("cli")
      suppressWarnings(with_progress({
        p <- progressr::progressor(steps = length(plots),
                                   message = 'compute stats')
        selpix <- mapply(FUN = get_samples_from_tiles,
                         plotID = names(plots), pix2sel = pix2sel,
                         MoreArgs = list(listfiles = listfiles,
                                         feat_list = feature_list,
                                         as.df = T,
                                         p = p),
                         SIMPLIFY = F)}))
    } else {
      message('compute stats')
      cl <- parallel::makeCluster(nbCPU)
      plan("cluster", workers = cl)
      selpix <- future.apply::future_mapply(FUN = get_samples_from_tiles,
                                            plotID = names(plots), pix2sel = pix2sel,
                                            MoreArgs = list(listfiles = listfiles,
                                                            feat_list = feature_list,
                                                            as.df = T,
                                                            p = p),
                                            future.seed = T, SIMPLIFY = F)
      parallel::stopCluster(cl)
      plan(sequential)
    }
    # compute IQR
    selpixAll <- do.call(what = 'rbind', selpix)
    iqr_si <- lapply(X = selpixAll, FUN = biodivMapR::IQR_outliers, weightIRQ = weightIRQ)
    
    ##############################################################################
    # produce mask for each tile
    if (nbCPU==1){
      handlers("cli")
      suppressWarnings(with_progress({
        
        p <- progressr::progressor(steps = length(plots),
                                   message = 'update mask')
        mask_path <- lapply(X = names(plots),
                            FUN = update_mask_from_tiles,
                            listfiles = listfiles,
                            iqr_si = iqr_si,
                            mask_dir = mask_dir,
                            p = p)}))
    } else {
      message('update mask')
      cl <- parallel::makeCluster(nbCPU)
      plan("cluster", workers = cl)
      mask_path <- future.apply::future_lapply(X = names(plots),
                                               FUN = update_mask_from_tiles,
                                               listfiles = listfiles,
                                               iqr_si = iqr_si,
                                               mask_dir = mask_dir,
                                               future.seed = TRUE)
      parallel::stopCluster(cl)
      plan(sequential)
    }
  }
  return(mask_path)
}
