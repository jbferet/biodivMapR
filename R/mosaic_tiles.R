#' compute spectral indices and composite for a given period of time
#'
#' @param pattern character. identify which files to mosaic
#' @param dir_path character. firectory where rasters to mosaic are stored
#' @param vrt_save character. where to save vrt
#' @param siteName character. name of the site
#' @param overwrite boolean
#'
#' @return mosaic_path
#' @importFrom terra vrt
#' @importFrom sf gdal_utils
#' @export
#'
mosaic_tiles <- function(pattern, dir_path, vrt_save, siteName = NULL,
                         overwrite = FALSE){
  # create vrt
  if (! is.null(siteName))
    siteName <- paste0('_', siteName, '_')
  listfiles <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  output_vrt_path <- file.path(getwd(), paste0(siteName, pattern,'_mosaic.vrt'))
  if (!file.exists(output_vrt_path))
    v <- terra::vrt(x = listfiles, filename = output_vrt_path)
  # create tiff from vrt
  mosaic_path <- file.path(dir_path, paste0(siteName, pattern,'_mosaic.tiff'))
  if (!file.exists(mosaic_path) | overwrite)
    sf::gdal_utils(util = 'translate', source = output_vrt_path,
                   destination = mosaic_path,
                   options = c("-co", "COMPRESS=LZW"))
  # create tiff from vrt
  dir.create(path = vrt_save, showWarnings = FALSE, recursive = TRUE)
  output_vrt_path2 <- file.path(vrt_save,
                                paste0(siteName, pattern,'_mosaic.vrt'))
  file.rename(from = output_vrt_path, to = output_vrt_path2)
  return(mosaic_path)
}
