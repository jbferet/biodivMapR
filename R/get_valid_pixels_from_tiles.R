#' get valid pixels from list of tiles
#'
#' @param feature_dir character. directory where features to be used by biodivMapR are stored.
#' @param plots list.
#' @param nbCPU numeric.
#' @param mask_dir character.
#'
#' @return nbPixValid
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply

#' @export

get_valid_pixels_from_tiles <- function(feature_dir, plots, nbCPU = 1,
                                        mask_dir = NULL){

  nbPixValid <- NULL
  # identify how many samples per tile should be extracted
  listfiles <- list.files(feature_dir, full.names = T)
  if (!is.null(mask_dir))
    listfiles <- list.files(mask_dir, full.names = T)
  ##############################################################################
  # get number of pixels per tile
  # if (nbCPU==1){
    handlers("cli")
    suppressWarnings(with_progress({
      p <- progressr::progressor(steps = length(plots),
                                 message = 'get valid pixels from tiles')
      nbPixValid <- lapply(X = names(plots), FUN = get_valid_pixels,
                           listfiles = listfiles, p = p)}))
  # } else {
  #   message('get valid pixels from tiles')
  #   cl <- parallel::makeCluster(nbCPU)
  #   plan("cluster", workers = cl)
  #   nbPixValid <- future.apply::future_lapply(X = names(plots),
  #                                             FUN = get_valid_pixels,
  #                                             listfiles = listfiles,
  #                                             future.seed = TRUE)
  #   parallel::stopCluster(cl)
  #   plan(sequential)
  # }
  names(nbPixValid) <- names(plots)
  return(nbPixValid)
}
