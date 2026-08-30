#' computes a mask for a collection of tiles based on interquartile range filtering
#' interquartile is computed for all data included in teh collection
#' samples out of a range defined by
#' quartile_1-weightIQR x (quartile_3 - quartile_1) and
#' quartile_3+weightIQR x (quartile_3-quartile_1)
#' are then considered as outliers
#'
#' @param feature_dir character. directory where all individual features are stored
#' @param feature_list character. list of features. files are expected to end with '_feature.'
#' @param mask_dir character.
#' @param plot_names character. plot names
#' @param weightIQR numeric.
#' @param filetype character.
#' @param nbCPU numeric.
#' @param nb_pixstats numeric.
#'
#' @return mask_path_list
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom parallel makeCluster stopCluster
#' @export
#'
compute_mask_iqr_tiles <- function(feature_dir, feature_list, mask_dir, plot_names,
                                   weightIQR = 4, filetype = 'GTiff', nbCPU = 1,
                                   nb_pixstats = 5e6){

  dir.create(path = mask_dir, showWarnings = F, recursive = T)
  maxCPU <- length(plot_names)
  if (nbCPU > maxCPU)
    nbCPU <-  maxCPU

  process_mask <- TRUE
  # first check missing masks
  mask_path <- file.path(mask_dir, paste0('mask_', plot_names, '_IQR.tiff'))
  tile_exists <- plot_names[which(file.exists(mask_path))]
  mask_missing0 <- plot_names[which(!file.exists(mask_path))]
  mask_missing <- paste0('_', plot_names[which(!file.exists(mask_path))], '_')
  if (length(mask_missing0)==0){
    process_mask <- FALSE
    tile_exists <- plot_names[which(file.exists(mask_path))]
  } else {
    # check if features also exist
    feature_list_extended <- paste0('_',feature_list, '.')
    if (filetype == 'GTiff')
      feature_list_extended <- paste0('_',feature_list, '.tiff')
    features_files <- lapply(X = feature_list_extended,
                             FUN = list.files, path = feature_dir)
    names(features_files) <- feature_list
    feat_exists <- list()
    # which features exist
    for (feat in feature_list)
      feat_exists[[feat]] <- unlist(lapply(X = mask_missing,
                                           FUN = grep, x = features_files[[feat]]))
    if (length(unlist(feat_exists))==0)
      process_mask <- FALSE
  }

  if (process_mask){
    # list files
    feature_list_extended <- paste0('_',feature_list, '.')
    if (filetype == 'GTiff')
      feature_list_extended <- paste0('_',feature_list, '.tiff')
    listfiles <- unlist(lapply(X = feature_list_extended,
                               FUN = list.files, path = feature_dir,
                               full.names = TRUE))
    ############################################################################
    # get number of pixels per tile
    # if (nbCPU==1){
    handlers("cli")
    suppressWarnings(with_progress({
      p <- progressr::progressor(steps = length(plot_names),
                                 message = 'get valid pixels from tiles')
      nb_pix_valid <- lapply(X = plot_names, FUN = get_valid_pixels,
                             listfiles = listfiles, p = p)}))
    # } else {
    #   message('get valid pixels from tiles')
    #   cl <- parallel::makeCluster(nbCPU)
    #   plan("cluster", workers = cl)
    #   nb_pix_valid <- future.apply::future_lapply(X = names(plots),
    #                                             FUN = get_valid_pixels,
    #                                             listfiles = listfiles,
    #                                             future.seed = TRUE)
    #   parallel::stopCluster(cl)
    #   plan(sequential)
    # }
    # get total number of pixels
    totalPixels <- sum(unlist(nb_pix_valid))
    ############################################################################
    # get statistics on nb_pixstats
    if (totalPixels<nb_pixstats)
      nb_pixstats <- totalPixels
    ratioStats <- nb_pixstats/totalPixels
    pix2sel <- lapply(X = nb_pix_valid,
                      FUN = function(x, ratio){as.numeric(round(x*ratio))},
                      ratio = ratioStats)
    if (nbCPU==1){
      handlers("cli")
      suppressWarnings(with_progress({
        p <- progressr::progressor(steps = length(plot_names),
                                   message = 'compute stats')
        selpix <- mapply(FUN = get_samples_from_tiles,
                         plotID = plot_names, pix2sel = pix2sel,
                         MoreArgs = list(listfiles = listfiles,
                                         feat_list = feature_list,
                                         as.df = TRUE,
                                         p = p),
                         SIMPLIFY = FALSE)}))
    } else {
      message('compute stats')
      cl <- parallel::makeCluster(nbCPU)
      with(plan("cluster", workers = cl), local = TRUE)
      selpix <- future.apply::future_mapply(FUN = get_samples_from_tiles,
                                            plotID = plot_names,
                                            pix2sel = pix2sel,
                                            MoreArgs = list(listfiles = listfiles,
                                                            feat_list = feature_list,
                                                            as.df = TRUE),
                                            future.seed = TRUE,
                                            future.chunk.size = NULL,
                                            future.scheduling = 1, SIMPLIFY = FALSE)
      parallel::stopCluster(cl)
      plan(sequential)
    }
    # compute IQR
    selpixAll <- do.call(what = 'rbind', selpix)
    iqr_si <- lapply(X = selpixAll,
                     FUN = iqr_outliers,
                     weightIQR = weightIQR)

    ##############################################################################
    # produce mask for each tile
    # if (nbCPU==1){
    handlers("cli")
    suppressWarnings(with_progress({

      p <- progressr::progressor(steps = length(plot_names),
                                 message = 'update mask')
      mask_path <- lapply(X = plot_names,
                          FUN = update_mask_from_tiles,
                          listfiles = listfiles,
                          iqr_si = iqr_si,
                          mask_dir = mask_dir,
                          p = p)}))
    # } else {
    #   message('update mask')
    #   cl <- parallel::makeCluster(nbCPU)
    #   plan("cluster", workers = cl)
    #   mask_path <- future.apply::future_lapply(X = names(plots),
    #                                            FUN = update_mask_from_tiles,
    #                                            listfiles = listfiles,
    #                                            iqr_si = iqr_si,
    #                                            mask_dir = mask_dir,
    #                                            future.seed = TRUE)
    #   parallel::stopCluster(cl)
    #   plan(sequential)
    # }
    notnullMask <- which(!unlist(lapply(X = mask_path, FUN = is.null)))
    tile_exists <- plot_names[notnullMask]
    mask_path <- mask_path[notnullMask]
  }
  mask_path_list <- list('mask_path' = mask_path,
                         'tile_exists' = tile_exists)
  return(mask_path_list)
}
