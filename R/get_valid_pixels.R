#' get valid pixels from a list of plots produced with 
#' preprocS2 function get_s2_tiling (notNA)
#' 
#' @param plotID list.
#' @param listfiles character.
#' @param p list.
#'
#' @return nbPixValid
#' @importFrom dplyr group_split %>%
#' @importFrom terra spatSample extract
#' @importFrom sf st_as_sf st_sf st_crs
#' @export
#'
get_valid_pixels <- function(plotID, listfiles, p = NULL){
  plotID <- paste0('_',plotID,'_')
  # tileSI <- listfiles[stringr::str_detect(string = listfiles, pattern = plotID)]
  tileSI <- listfiles[grepl(x = listfiles, pattern = plotID)]
  # get statistics on data availability
  nbPixValid <- 0
  if (length(tileSI) > 0){
    if (file.exists(tileSI[[1]])){
      rastID <- terra::rast(tileSI[[1]])
      nbPixValid <- terra::global(x = rastID, fun = 'notNA')
    }
  }
  if (!is.null(p)) p()
  return(nbPixValid)
}
