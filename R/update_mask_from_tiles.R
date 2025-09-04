#' update mask based on IQR for a series of rasters
#' produced with preprocS2 function get_s2_tiling
#'
#' @param plotID character.
#' @param listfiles character.
#' @param iqr_si character.
#' @param mask_dir list.
#' @param p numeric.
#'
#' @return filename
#' @importFrom terra values sources rast writeRaster
#' @export
#'
update_mask_from_tiles <- function(plotID, listfiles, iqr_si, mask_dir,
                                   p = NULL){
  plotID <- paste0('_',plotID,'_')
  # tileSI <- listfiles[stringr::str_detect(string = listfiles, pattern = plotID)]
  tileSI <- listfiles[grepl(x = basename(listfiles), pattern = plotID)]
  # get statistics on data availability
  namefeatures <- names(iqr_si)
  filename <- NULL
  if (length(tileSI) > 0){
    mask <- 1+0*terra::rast(tileSI[[1]])
    names(mask) <- mask
    # mask[is.na(mask)] <- 0
    rastSI <- terra::rast(tileSI)
    for (feat in namefeatures){
      # whichfeat <- which(stringr::str_detect(basename(terra::sources(rastSI)), feat) )
      whichfeat <- which(grepl(x = basename(terra::sources(rastSI)),
                               pattern = feat))
      names(rastSI)[whichfeat] <- feat
    }
    for (feat in namefeatures){
      elim <- which(terra::values(rastSI[[feat]])<iqr_si[[feat]][1] |
                      terra::values(rastSI[[feat]])>iqr_si[[feat]][2])
      if (length(elim)>0)
        mask[elim] <- NA
    }
    filename <- file.path(mask_dir, paste0('mask', plotID, 'IQR.tiff'))
    terra::writeRaster(x = mask, filename = filename, filetype = 'GTiff',
                       overwrite = TRUE)
  }
  if (!is.null(p))
    p()
  return(filename)
}
