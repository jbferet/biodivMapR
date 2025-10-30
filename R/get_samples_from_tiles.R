#' get samples for alpha diversity mapping
#'
#' @param plotID list.
#' @param pix2sel numeric.
#' @param listfiles character.
#' @param feat_list list.
#' @param as.df boolean.
#' @param as.points boolean.
#' @param xy boolean.
#' @param method character. method for terra::spatSample ('random' or 'regular')
#' @param p list.
#'
#' @return selpix
#' @importFrom dplyr group_split %>%
#' @importFrom terra spatSample extract
#' @importFrom sf st_as_sf st_sf st_crs
#' @export
#'
get_samples_from_tiles <- function(plotID, pix2sel, listfiles, feat_list,
                                   as.df = FALSE, as.points = FALSE, xy = FALSE,
                                   method = 'regular', p = NULL){
  plotID <- paste0('_',plotID,'_')
  tileSI <- listfiles[grepl(x = basename(listfiles), pattern = plotID)]
  # get statistics on data availability
  selpix <- NULL
  if (length(tileSI) > 0 & pix2sel >0){
    if (all(file.exists(tileSI))){
      rastID <- terra::rast(tileSI)
      for (feat in feat_list){
        feat2 <- feat
        if (!feat == 'mask')
          feat2 <- paste0('_',feat, '.')
        whichfeat <- which(grepl(x = basename(terra::sources(rastID)),
                                 pattern = feat2))
        names(rastID)[whichfeat] <- feat
      }
      selpix <- terra::spatSample(x = rastID, size = as.numeric(pix2sel),
                                  method = method, na.rm = TRUE, as.df = as.df,
                                  as.points = as.points, xy = xy, warn = FALSE)
    }
  }
  if (!is.null(p))
    p()
  return(selpix)
}
