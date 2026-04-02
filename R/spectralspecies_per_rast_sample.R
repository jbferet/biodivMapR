#' get spectral species corresponding to raster sample extracted with extract
#
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param rast_sample dataframe.
#' @param selected_bands numeric. bands selected from input data
#
#' @return list
#' @importFrom tidyr nest
#' @export
#'
spectralspecies_per_rast_sample <- function(Kmeans_info, rast_sample,
                                            selected_bands = NULL){

  SSValid <- NULL
  nb_pix_per_plot <- data.frame(table(rast_sample$ID))
  # only get common plots between nb_pix_per_plot and nb_pix_per_plot_init
  if (length(nb_pix_per_plot$Freq)>0){
    if (is.null(selected_bands))
      selected_bands <- seq_len(dim(rast_sample)[2]-1)
    rast_sample_noID <- rast_sample
    rast_sample_noID$ID <- NULL
    SSValid <- get_spectralSpecies(inputdata = rast_sample_noID,
                                   Kmeans_info = Kmeans_info,
                                   selected_bands = selected_bands)
    SSValid$win_ID <- rast_sample$ID
  }
  return(list('SSValid' = SSValid))
}
