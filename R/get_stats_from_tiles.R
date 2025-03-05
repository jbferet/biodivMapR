#' get samples from tiles produced with preprocS2 function get_s2_tiling
#' in order to produce stats
#'
#' @param plotID list.
#' @param pix2sel numeric.
#' @param listfiles character.
#' @param SI_list list.
#' @param p list.
#'
#' @return selpix
#' @importFrom dplyr group_split %>%
#' @importFrom terra spatSample extract
#' @importFrom sf st_as_sf st_sf st_crs
#' @export
#'
get_stats_from_tiles <- function(plotID, pix2sel, listfiles, SI_list, p = NULL){
  plotID <- paste0('_',plotID,'_')
  # tileSI <- listfiles[stringr::str_detect(string = listfiles, pattern = plotID)]
  tileSI <- listfiles[grepl(x = listfiles, pattern = plotID)]
  # get statistics on data availability
  selpix <- NULL
  if (length(tileSI) > 0 & pix2sel >0){
    if (all(file.exists(tileSI))){
      rastID <- terra::rast(tileSI)
      for (feat in SI_list){
        # whichfeat <- which(stringr::str_detect(basename(terra::sources(rastID)), feat) )
        whichfeat <- which(grepl(x = basename(terra::sources(rastID)),
                                 pattern = feat))
        names(rastID)[whichfeat] <- feat
      }
      selpix <- spatSample(x = rastID, size = as.numeric(pix2sel), method = "random",
                           na.rm = T, as.df = T, warn = F)
    }
  }
  if (!is.null(p)) p()
  return(selpix)
}
