#'save diversity maps as raster data
#'
#' @param input_raster_path list. list of paths corresponding to input rasters
#' @param ab_div_metrics list produced from get_raster_diversity
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson, hill
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

save_diversity_maps_tile <- function(input_raster_path,
                                     ab_div_metrics,
                                     alphametrics = 'shannon',
                                     Hill_order = 1,
                                     FDmetric = NULL,
                                     input_rast,
                                     output_dir,
                                     output_raster_name = NULL,
                                     window_size,
                                     filetype = 'GTiff'){
  # save alpha diversity indices

  # rast_out <- terra::rast(input_raster_path[[1]])
  # for (idx in alphametrics){
  #   for (crit in c('mean', 'sd')){
  #     elem <- paste0(idx,'_',crit)
  #     names(rast_out) <- paste0(idx,'_',crit)
  #     terra::values(rast_out) <- res_shapeChunk[[idx]][[crit]]
  #     # save file
  #   }
  # }

  for (idx in alphametrics) {
    idx2 <- idx
    if (idx == 'hill') idx2 <- paste0(idx, '_', Hill_order)
    # Mean value
    # produce a template
    template_rast <- terra::rast(input_raster_path[[1]])
    template_newres <- terra::aggregate(x = template_rast, fact = window_size)
    terra::values(template_newres) <- c(ab_div_metrics[[idx]]$mean)
    names(template_newres[[1]]) <- paste('mean', idx)
    # define output raster name
    if (is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir, paste0(idx2, '_mean'))
    if (!is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir, paste0(output_raster_name[[idx2]], '_mean'))

    if (filetype%in%c('GTiff', 'COG')) output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_newres, filename = output_raster,
                       filetype = filetype, overwrite = T, gdal=c("COMPRESS=LZW"))
    # SD value
    # produce a template
    terra::values(template_newres[[1]]) <- c(ab_div_metrics[[idx]]$sd)
    names(template_newres[[1]]) <- paste('sd', idx)
    # define output raster name
    if (is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir, paste0(idx2, '_sd'))
    if (!is.null(output_raster_name[[idx2]]))
      output_raster <- file.path(output_dir, paste0(output_raster_name[[idx2]], '_sd'))
    if (filetype%in%c('GTiff', 'COG')) output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_newres, filename = output_raster,
                       filetype = filetype, overwrite = T, gdal=c("COMPRESS=LZW"))
  }

  # save functional diversity indices
  for (idx in FDmetric){
    # produce a template
    template_rast <- terra::rast(input_raster_path[[1]])
    template_newres <- terra::aggregate(x = template_rast, fact = window_size)
    terra::values(template_newres) <- c(ab_div_metrics[[idx]])
    names(template_newres[[1]]) <- idx
    # define output raster name
    if (is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir, idx)
    if (!is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir, output_raster_name[[idx]])
    if (filetype%in%c('GTiff', 'COG')) output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = T, gdal=c("COMPRESS=LZW"))
  }

  # save beta diversity indices
  # get PCoA and corresponding dimensions
  PCoA <- ab_div_metrics$PCoA_BC
  dimPCoA <- length(PCoA)
  # produce a template with N dimensions
  template_rast <- terra::rast(input_raster_path[[1]])
  template_newres <- terra::aggregate(x = template_rast, fact = window_size)
  dim(template_newres)[3] <- dimPCoA
  for (nbPC in seq_len(dimPCoA))
    terra::values(template_newres[[nbPC]]) <- c(PCoA[[nbPC]])
  names(template_newres) <- paste0('PCoA#',seq(1,dimPCoA))

  if (is.null(output_raster_name[['beta']]))
    output_raster <- file.path(output_dir, 'beta')
  if (!is.null(output_raster_name[[idx2]]))
    output_raster <- file.path(output_dir, output_raster_name[['beta']])
  if (filetype%in%c('GTiff', 'COG')) output_raster <- paste0(output_raster, '.tiff')
  terra::writeRaster(x = template_newres, filename = output_raster,
                     filetype = filetype, overwrite = T, gdal=c("COMPRESS=LZW"))
  return(invisible())
}
