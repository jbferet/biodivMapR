#'save diversity maps as raster data
#'
#' @param ab_div_metrics list produced from get_raster_diversity
#' @param alphametrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#' @param FDmetric character. list of functional metrics
#' @param input_rast list. path for input rasters
#' @param output_dir character.
#' @param output_raster_name character.
#' @param window_size numeric. window size for square plots
#' @param filetype character. GDAL driver
#'
#' @return blk
#' @importFrom terra aggregate writeRaster values names
#' @export

save_diversity_maps <- function(ab_div_metrics,
                                alphametrics = 'shannon',
                                Hill_order = 1,
                                FDmetric = NULL,
                                input_rast,
                                output_dir,
                                output_raster_name = NULL,
                                window_size,
                                filetype = 'GTiff'){

  diversity_maps <- list()
  # save alpha diversity indices
  for (idx in alphametrics) {
    idx2 <- idx
    if (idx == 'hill')
      idx2 <- paste0(idx, '_', Hill_order)
    # Mean value
    # produce a template
    template_rast <- terra::aggregate(input_rast[[1]], fact = window_size)
    mat_idx_mean <- array(do.call(rbind,lapply(lapply(ab_div_metrics,
                                                      '[[', idx),
                                               '[[','mean')),
                          dim = dim(template_rast))
    terra::values(template_rast[[1]]) <- mat_idx_mean
    names(template_rast[[1]]) <- paste('mean', idx)
    # define output raster name
    if (is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir, paste0(idx2, '_mean'))
    if (!is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir,
                                 paste0(output_raster_name[[idx2]], '_mean'))

    if (filetype%in%c('GTiff', 'COG'))
      output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW"))

    diversity_maps[[paste0(idx2, '_mean')]] <- output_raster

    # SD value
    # produce a template
    template_rast <- terra::aggregate(input_rast[[1]], fact = window_size)
    mat_idx_sd <- array(do.call(rbind,lapply(lapply(ab_div_metrics, '[[', idx),
                                             '[[','sd')),
                        dim = dim(template_rast))
    terra::values(template_rast[[1]]) <- mat_idx_sd
    names(template_rast[[1]]) <- paste('sd', idx)
    # define output raster name
    if (is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir, paste0(idx2, '_sd'))
    if (!is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir, paste0(output_raster_name[[idx2]],
                                                    '_sd'))
    if (filetype%in%c('GTiff', 'COG'))
      output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW"))
    diversity_maps[[paste0(idx2, '_sd')]] <- output_raster
  }

  # save functional diversity indices
  for (idx in FDmetric){
    # produce a template
    template_rast <- terra::aggregate(input_rast[[1]], fact = window_size)
    mat_idx <- array(do.call(rbind,lapply(ab_div_metrics, '[[', idx)),
                     dim = dim(template_rast))
    terra::values(template_rast[[1]]) <- mat_idx
    names(template_rast[[1]]) <- idx
    # define output raster name
    if (is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir, idx)
    if (!is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir, output_raster_name[[idx]])
    if (filetype%in%c('GTiff', 'COG'))
      output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW"))
    diversity_maps[[idx]] <- output_raster
  }

  # save beta diversity indices
  # get PCoA and corresponding dimensions
  PCoA <- lapply(ab_div_metrics, '[[', 'PCoA_BC')
  dimPCoA <- dim(PCoA[[1]])[3]
  # produce a template with N dimensions
  template_rast <- terra::aggregate(input_rast[[1]], fact = window_size)
  dim(template_rast)[3] <- dimPCoA
  for (nbPC in seq_len(dimPCoA)){
    pctmp <- list()
    for (nbPieces in seq_len(length(PCoA)))
      pctmp[[nbPieces]] <- PCoA[[nbPieces]][,,nbPC]
    terra::values(template_rast[[nbPC]]) <- do.call(rbind,pctmp)
  }
  names(template_rast) <- paste0('PCoA#',seq(1,dimPCoA))

  if (is.null(output_raster_name[['beta']]))
    output_raster <- file.path(output_dir, 'beta')
  if (!is.null(output_raster_name[[idx2]]))
    output_raster <- file.path(output_dir, output_raster_name[['beta']])
  if (filetype%in%c('GTiff', 'COG'))
    output_raster <- paste0(output_raster, '.tiff')
  terra::writeRaster(x = template_rast, filename = output_raster,
                     filetype = filetype, overwrite = TRUE,
                     gdal = c("COMPRESS=LZW"))
  diversity_maps[['beta']] <- output_raster

  return(diversity_maps)
}
