#' mosaic tiles defined in vrt and compress if possible
#'
#' @param vrt_path character. path for vrt file
#' @param mosaic_path character. path for mosaic file
#'
#' @return mosaic_path
#' @importFrom sf gdal_utils
#' @export
#'
mosaic_from_vrt <- function(vrt_path, mosaic_path){
  result <- try({
    sf::gdal_utils(util = 'translate', source = vrt_path,
                   destination = mosaic_path,
                   options = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"))
    # co = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"))
  }, silent = TRUE)
  if ( "try-error" %in% class(result) ) {
    sf::gdal_utils(util = 'translate', source = vrt_path,
                   destination = mosaic_path)
  }
}
