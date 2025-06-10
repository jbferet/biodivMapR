#' apply biodivMapR to a set of plots identified by a field '_id_'
#' produced with preprocS2 function get_s2_tiling
#'
#' @param id character. ID for plot
#' @param feature_dir directory where spectral indices are for each plot
#' @param mask_dir character.
#' @param list_features character.
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param Beta_info list. BC dissimilarity & associated beta metrics
#' @param output_dir path where to save outputs
#' @param selected_bands numeric. bands selected from input data
#' @param window_size numeric. window size for square plots
#' @param alphametrics list. alpha diversity metrics
#' @param Hill_order numeric. Hill order
#' @param FDmetric character. list of functional metrics
#' @param pcelim numeric. min proportion of pix to consider spectral species
#' @param nbCPU numeric. Number of CPUs available
#' @param maxRows numeric. maximum number or rows to process
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#' @param filetype character. gdal driver for output raster
#' @param moving_window boolean. should process be moving window (much longer)
#' @param p list. progressbar
#'
#' @return none
#' @export

run_biodivMapR_plot <- function(id, feature_dir, mask_dir = NULL,
                                list_features, Kmeans_info, Beta_info,
                                output_dir, selected_bands = NULL, window_size,
                                alphametrics = 'shannon', Hill_order = 1,
                                FDmetric = NULL, pcelim = 0.02,
                                maxRows = NULL, nbCPU = 1, min_sun = 0.25,
                                filetype = 'GTiff', moving_window = FALSE,
                                p = NULL){

  betanames <- paste0('beta_',id)
  alphanames <- paste0(alphametrics,'_',id)
  functionalname <- NULL
  if (!is.null(FDmetric))
    functionalname <- paste0(FDmetric,'_',id)
  alphanames_mean <- paste0(alphanames,'_mean')
  output_raster_name <- as.list(c(betanames, alphanames, functionalname))
  output_raster_name_mean <- as.list(c(betanames, alphanames_mean, functionalname))
  names(output_raster_name) <- c('beta', alphametrics, FDmetric)
  if (FALSE %in% file.exists(file.path(output_dir,
                                       paste0(output_raster_name_mean,'.tiff')))){
    list_feat <- list.files(path = feature_dir, pattern = paste0('_',id,'_'))
    list_feat <- unique(gsub(pattern = '.aux.xml', replacement = '',
                             x = list_feat))
    if (!is.null(mask_dir)){
      input_mask_path <- list.files(path = mask_dir,
                                    pattern = paste0('_',id,'_'),
                                    full.names = TRUE)
      input_mask_path <- unique(gsub(pattern = '.aux.xml', replacement = '',
                                     x = input_mask_path))
    } else {
      input_mask_path <- NULL
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
      if (!FALSE %in% file.exists(unlist(input_raster_path)))
        run_biodivMapR(input_raster_path = input_raster_path,
                       input_mask_path = input_mask_path,
                       Kmeans_info = Kmeans_info, Beta_info = Beta_info,
                       output_dir = output_dir,
                       output_raster_name = output_raster_name,
                       selected_bands = selected_bands,
                       window_size = window_size,
                       alphametrics = alphametrics, Hill_order = Hill_order,
                       FDmetric = FDmetric, pcelim = pcelim,
                       maxRows = maxRows, nbCPU = nbCPU, min_sun = min_sun,
                       filetype = filetype, moving_window = moving_window)
    }
  }
  if (!is.null(p))
    p()
  return()
}
