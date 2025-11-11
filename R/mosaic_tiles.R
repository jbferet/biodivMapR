#' compute spectral indices and composite for a given period of time
#'
#' @param pattern character. identify which files to mosaic
#' @param dir_path character. firectory where rasters to mosaic are stored
#' @param vrt_save character. where to save vrt
#' @param site_name character. name of the site
#' @param overwrite boolean
#'
#' @return mosaic_path
#' @importFrom terra vrt
#' @importFrom sf gdal_utils
#' @export
#'
mosaic_tiles <- function(pattern, dir_path, vrt_save, site_name = NULL,
                         overwrite = FALSE){
  # create vrt
  if (! is.null(site_name))
    site_name <- paste0(site_name, '_')
    # site_name <- paste0('_', site_name, '_')
  listfiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  output_vrt_path <- file.path(getwd(),
                               paste0(site_name, pattern,'_mosaic.vrt'))
  # if (!file.exists(output_vrt_path))
  v <- terra::vrt(x = listfiles, filename = output_vrt_path, overwrite = TRUE)

  # create tiff from vrt
  mosaic_path <- file.path(dir_path, paste0(site_name, pattern,'_mosaic.tiff'))
  if (!file.exists(mosaic_path) | overwrite){
    message(paste('write image for diversity metric', pattern))
    result <- try({
      sf::gdal_utils(util = 'translate', source = output_vrt_path,
                     destination = mosaic_path,
                     options = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"))
      # co = c("COMPRESS=LZW", "BIGTIFF=IF_SAFER"))
    }, silent = TRUE)
    if ( "try-error" %in% class(result) ) {
      sf::gdal_utils(util = 'translate', source = output_vrt_path,
                     destination = mosaic_path)
    }
  }
  # delete vrt
  file.remove(output_vrt_path)
  # delete individual files
  file.remove(listfiles)
  return(mosaic_path)
}
