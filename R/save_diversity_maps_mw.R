#'save diversity maps as raster data once moving window process done
#'
#' @param input_raster_path list. list of paths corresponding to input rasters
#' @param ab_div_metrics list produced from get_raster_diversity
#' @param alpha_metrics list. alpha diversity metrics: richness, shannon, simpson, hill
#' @param Hill_order numeric. Hill order
#' @param beta_metrics boolean. set TRUE to compute beta diversity
#' @param fd_metrics character. list of functional metrics
#' @param input_rast list. path for input rasters
#' @param output_dir character.
#' @param output_raster_name character.
#' @param window_size numeric. window size for square plots
#' @param filetype character. GDAL driver
#'
#' @return blk
#' @importFrom terra aggregate writeRaster values names
#' @export

save_diversity_maps_mw <- function(input_raster_path,
                                   ab_div_metrics,
                                   alpha_metrics = 'shannon',
                                   Hill_order = 1,
                                   beta_metrics = TRUE,
                                   fd_metrics = NULL,
                                   input_rast,
                                   output_dir,
                                   output_raster_name = NULL,
                                   window_size,
                                   filetype = 'GTiff'){
  # save alpha diversity indices

  # rast_out <- terra::rast(input_raster_path[[1]])
  # for (idx in alpha_metrics){
  #   for (crit in c('mean', 'sd')){
  #     elem <- paste0(idx,'_',crit)
  #     names(rast_out) <- paste0(idx,'_',crit)
  #     terra::values(rast_out) <- res_shapeChunk[[idx]][[crit]]
  #     # save file
  #   }
  # }

  diversity_maps <- list()
  for (idx in alpha_metrics) {
    idx2 <- idx
    if (idx == 'hill')
      idx2 <- paste0(idx, '_', Hill_order)
    # Mean value
    # produce a template
    template_rast <- terra::rast(input_raster_path[[1]])
    terra::values(template_rast[[1]]) <- ab_div_metrics[[idx]]$mean
    names(template_rast[[1]]) <- paste('mean', idx)
    # crop borders to avoid artifacts
    template_rast <- crop_border_tile(rast_obj = template_rast,
                                      window_size = window_size)
    # define output raster name
    if (is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir, paste0(idx2, '_mean'))
    if (!is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir,
                                 paste0(output_raster_name[[idx]], '_mean'))
    if (filetype%in%c('GTiff', 'COG'))
      output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW"))
    diversity_maps[[paste0(idx2, '_mean')]] <- output_raster
    # SD value
    # produce a template
    template_rast <- terra::rast(input_raster_path[[1]])
    terra::values(template_rast[[1]]) <- ab_div_metrics[[idx]]$sd
    names(template_rast[[1]]) <- paste('sd', idx)
    # crop borders to avoid artifacts
    template_rast <- crop_border_tile(rast_obj = template_rast,
                                      window_size = window_size)
    # define output raster name
    if (is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir, paste0(idx2, '_sd'))
    if (!is.null(output_raster_name[[idx]]))
      output_raster <- file.path(output_dir,
                                 paste0(output_raster_name[[idx]], '_sd'))
    if (filetype%in%c('GTiff', 'COG'))
      output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW"))
    diversity_maps[[paste0(idx2, '_sd')]] <- output_raster
  }

  # save functional diversity indices
  for (idx in fd_metrics){
    # produce a template
    template_rast <- terra::rast(input_raster_path[[1]])
    terra::values(template_rast[[1]]) <- ab_div_metrics[[idx]]
    names(template_rast[[1]]) <- idx
    # crop borders to avoid artifacts
    template_rast <- crop_border_tile(rast_obj = template_rast,
                                      window_size = window_size)
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
  if (beta_metrics){
    # get PCoA and corresponding dimensions
    PCoA <- ab_div_metrics$PCoA_BC
    dimPCoA <- length(PCoA)
    # produce a template with N dimensions
    template_rast <- terra::rast(input_raster_path[[1]])
    dim(template_rast)[3] <- dimPCoA
    for (nbPC in seq_len(dimPCoA))
      terra::values(template_rast[[nbPC]]) <- PCoA[[nbPC]]
    names(template_rast) <- paste0('PCoA#',seq(1,dimPCoA))
    # crop borders to avoid artifacts
    template_rast <- crop_border_tile(rast_obj = template_rast,
                                      window_size = window_size)

    if (is.null(output_raster_name[['beta']]))
      output_raster <- file.path(output_dir, 'beta')
    if (!is.null(output_raster_name[['beta']]))
      output_raster <- file.path(output_dir, output_raster_name[['beta']])
    if (filetype%in%c('GTiff', 'COG'))
      output_raster <- paste0(output_raster, '.tiff')
    terra::writeRaster(x = template_rast, filename = output_raster,
                       filetype = filetype, overwrite = TRUE,
                       gdal = c("COMPRESS=LZW"))
    diversity_maps[['beta']] <- output_raster
  }
  return(diversity_maps)
}
