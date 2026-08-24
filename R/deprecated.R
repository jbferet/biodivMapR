#' Deprecated functions in biodivMapR
#'
#' These functions are provided for compatibility with older versions of biodivMapR but are deprecated.
#'
#' @param input_raster_path See the new functions that replace the old ones
#' @param input_mask_path See the new functions that replace the old ones
#' @param output_dir See the new functions that replace the old ones
#' @param window_size See the new functions that replace the old ones
#' @param selected_bands See the new functions that replace the old ones
#' @param Kmeans_info_save See the new functions that replace the old ones
#' @param Kmeans_info_read See the new functions that replace the old ones
#' @param Beta_info_save See the new functions that replace the old ones
#' @param Beta_info_read See the new functions that replace the old ones
#' @param nbCPU See the new functions that replace the old ones
#' @param options See the new functions that replace the old ones
#' @param feature_dir See the new functions that replace the old ones
#' @param list_features See the new functions that replace the old ones
#' @param mask_dir See the new functions that replace the old ones
#' @param plots See the new functions that replace the old ones
#' @param site_name See the new functions that replace the old ones
#'
#' @rdname deprecated
#' @name deprecated
#'
NULL

#' @export
#' @rdname deprecated
biodivMapR_full <- function(
    input_raster_path, input_mask_path = NULL, output_dir, window_size,
    selected_bands = NULL, Kmeans_info_save = NULL, Kmeans_info_read = NULL,
    Beta_info_save = NULL, Beta_info_read = NULL,  nbCPU = 1, options = NULL){

  .Deprecated("biodivMapR")
  return(biodivMapR(input_raster_path = input_raster_path,
                    input_mask_path = input_mask_path, output_dir = output_dir,
                    window_size = window_size, selected_bands = selected_bands,
                    Kmeans_info_save = Kmeans_info_save, Kmeans_info_read = Kmeans_info_read,
                    Beta_info_save = Beta_info_save, Beta_info_read = Beta_info_read,
                    nbCPU = nbCPU, options = options))
}


#' @export
#' @rdname deprecated
biodivMapR_full_tiles <- function(
    feature_dir, list_features, mask_dir = NULL, output_dir, window_size, plots,
    nbCPU = 1, site_name = NULL, options = NULL){

  .Deprecated("biodivMapR_tiling")
  return(biodivMapR_tiling(feature_dir = feature_dir,
                           list_features = list_features,
                           mask_dir = mask_dir,
                           output_dir = output_dir,
                           window_size = window_size,
                           plots = plots,
                           nbCPU = nbCPU,
                           site_name = site_name,
                           options = options))
}


