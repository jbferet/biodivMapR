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
                                   as.df = F, as.points = F, xy = F,
                                   method = 'regular', p = NULL){
  plotID <- paste0('_',plotID,'_')
  # tileSI <- listfiles[stringr::str_detect(string = listfiles, pattern = plotID)]
  tileSI <- listfiles[grepl(x = listfiles, pattern = plotID)]
  # get statistics on data availability
  selpix <- NULL
  if (length(tileSI) > 0 & pix2sel >0){
    if (all(file.exists(tileSI))){
      rastID <- terra::rast(tileSI)
      for (feat in feat_list){
        # whichfeat <- which(stringr::str_detect(basename(terra::sources(rastID)), feat))
        whichfeat <- which(grepl(x = basename(terra::sources(rastID)),
                                 pattern = feat))
        names(rastID)[whichfeat] <- feat
      }
      selpix <- spatSample(x = rastID, size = as.numeric(pix2sel), method = method,
                           na.rm = T, as.df = as.df, as.points = as.points,
                           xy = xy, warn = F)
    }
  }
  if (!is.null(p)) p()
  return(selpix)
}
