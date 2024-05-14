#'save diversity maps as raster data
#'
#' @param ab_div_metrics list produced from get_raster_diversity
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param input_rast list. path for input rasters
#' @param output_dir character.
#' @param window_size numeric. window size for square plots
#' @param filetype character. GDAL driver
#'
#' @return blk
#' @importFrom terra aggregate writeRaster values names
#' @export

save_diversity_maps <- function(ab_div_metrics,
                                alphametrics = 'shannon',
                                input_rast,
                                output_dir,
                                window_size,
                                filetype = 'COG'){
  # save alpha diversity indices
  for (idx in alphametrics) {
    # Mean value
    # produce a template
    template_rast <- terra::aggregate(input_rast[[1]], fact = window_size)
    mat_idx_mean <- array(do.call(rbind,lapply(lapply(ab_div_metrics, '[[', idx),
                                               '[[','mean')),
                          dim = dim(template_rast))
    terra::values(template_rast[[1]]) <- mat_idx_mean
    names(template_rast[[1]]) <- paste('mean', idx)
    # define output raster name
    output_raster <- file.path(output_dir, paste0(idx, '_mean'))
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = T)
    # SD value
    # produce a template
    template_rast <- terra::aggregate(input_rast[[1]], fact = window_size)
    mat_idx_sd <- array(do.call(rbind,lapply(lapply(ab_div_metrics, '[[', idx),
                                             '[[','sd')),
                        dim = dim(template_rast))
    terra::values(template_rast[[1]]) <- mat_idx_sd
    names(template_rast[[1]]) <- paste('sd', idx)
    # define output raster name
    output_raster <- file.path(output_dir, paste0(idx, '_sd'))
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = T)
  }

  # save beta diversity indices
  # get PCoA and corresponding dimensions
  PCoA <- lapply(ab_div_metrics, '[[', 'PCoA_BC')
  dimPCoA <- dim(PCoA[[1]])[3]
  # produce a template with N dimensions
  template_rast <- terra::aggregate(input_rast[[1]], fact = window_size)
  dim(template_rast)[3] <- dimPCoA
  for (nbPC in 1:dimPCoA){
    pctmp <- list()
    for (nbPieces in 1:length(PCoA)) pctmp[[nbPieces]] <- PCoA[[nbPieces]][,,nbPC]
    terra::values(template_rast[[nbPC]]) <- do.call(rbind,pctmp)
  }
  names(template_rast) <- paste0('PCoA#',seq(1,dimPCoA))
  output_raster <- file.path(output_dir, 'beta')
  terra::writeRaster(x = template_rast, filename = output_raster,
                     filetype = filetype, overwrite = T)
  return(invisible())
}
