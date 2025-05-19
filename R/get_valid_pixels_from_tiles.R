#' get valid pixels from list of tiles
#'
#' @param feature_dir character. directory where features to be used by biodivMapR are stored.
#' @param plots list.
#' @param nbCPU numeric.
#' @param mask_dir character.
#'
#' @return nb_pix_valid
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply

#' @export

get_valid_pixels_from_tiles <- function(feature_dir, plots, nbCPU = 1,
                                        mask_dir = NULL){

  nb_pix_valid <- NULL
  # identify how many samples per tile should be extracted
  listfiles <- list.files(feature_dir, full.names = TRUE)
  if (!is.null(mask_dir))
    listfiles <- list.files(mask_dir, full.names = TRUE)
  # get number of pixels per tile
  handlers("cli")
  suppressWarnings(with_progress({
    p <- progressr::progressor(steps = length(plots),
                               message = 'get valid pixels from tiles')
    nb_pix_valid <- lapply(X = names(plots), FUN = get_valid_pixels,
                           listfiles = listfiles, p = p)}))
  names(nb_pix_valid) <- names(plots)
  return(nb_pix_valid)
}
