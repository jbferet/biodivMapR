#' apply biodivMapR to a set of plots identified by a field '_id_'
#' produced with preprocS2 function get_s2_tiling
#'
#' @param id character. ID for plot
#' @param feature_dir directory where spectral indices are for each plot
#' @param mask_dir character.
#' @param list_features character.
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param output_dir path where to save outputs
#' @param selected_bands numeric. bands selected from input data
#' @param filetype character. gdal driver for output raster
#' @param overwrite boolean.
#' @param p list. progressbar
#'
#' @return none
#' @export

apply_spectral_species_plot <- function(id, feature_dir, mask_dir = NULL,
                                        list_features, Kmeans_info, output_dir,
                                        selected_bands = NULL, filetype = 'GTiff',
                                        overwrite = TRUE, p = NULL){

  output_raster_name <- as.list(paste0('spectral_species_',id))
  names(output_raster_name) <- 'spectral_species'
  output_raster_full_name <- file.path(output_dir, paste0(output_raster_name,'.tiff'))
  if (FALSE %in% file.exists(output_raster_full_name) | overwrite){
    list_feat <- list.files(path = feature_dir, pattern = paste0('_',id,'_'))
    list_feat <- clean_tiff_list(list_feat)
    input_mask_path <- NULL
    if (!is.null(mask_dir)){
      input_mask_path <- list.files(path = mask_dir,
                                    pattern = paste0('_',id,'_'),
                                    full.names = TRUE)
      input_mask_path <- clean_tiff_list(input_mask_path)
    }
    selfeat <- c()
    for (feat in list_features)
      selfeat <- c(selfeat, grep(x = list_feat,
                                 pattern = paste0('_',feat,'.tiff')))
    if (length(selfeat) == length(list_features)){
      feature_files <- list_feat[selfeat]
      names(feature_files) <- list_features
      input_raster_path <- list()
      for (feat in list_features)
        input_raster_path[[feat]] <- file.path(feature_dir, feature_files[feat])
      # prepare to read input raster data
      r_in <- list()
      if (is.null(names(input_raster_path)))
        names(input_raster_path) <- seq_len(length(input_raster_path))
      for (fid in names(input_raster_path))
        r_in[[fid]] <- terra::rast(input_raster_path[[fid]])
      # if a mask file is provided
      if (!is.null(input_mask_path)) {
        r_in[['mask']] <- terra::rast(input_mask_path)
        names(r_in[['mask']]) <- 'mask'
        input_raster_path[['mask']] <- input_mask_path
      }
      rast_in <- lapply(X = input_raster_path, FUN = terra::rast)
      rast_in <- terra::rast(rast_in)
      selected_bands <- seq_len(dim(rast_in)[3]-1)
      if (is.null(input_mask_path))
        selected_bands <- seq_len(dim(rast_in)[3])
      inputdata <- as.data.frame(terra::values(rast_in))
      sel <- which(inputdata$mask==1)
      inputdata <- inputdata[sel,]
      ss_rast <- NA*rast_in[[1]]
      names(ss_rast) <- 'spectral_species'
      if (length(inputdata)>0){
        ss_tile <- get_spectralSpecies(inputdata = inputdata,
                                       Kmeans_info = Kmeans_info,
                                       selected_bands = selected_bands)
        ss_rast[sel] <- ss_tile[[1]]
      }
      terra::writeRaster(x = ss_rast, filename = output_raster_full_name,
                         filetype = filetype, overwrite = TRUE,
                         gdal = c("COMPRESS=LZW"))
    }
  }
  if (!is.null(p))
    p()
  return(output_raster_full_name)
}

