#' get spectral species corresponding to raster sample extracted with extract
#
#' @param input_rast SpatRaster.
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param selected_bands numeric. bands selected from input data
#' @param fd_metrics boolean. should functional diversity be computed as well?
#' @param input_mask SpatRaster.
#' @param rast_sample dataframe.
#' @param AttributeTable dataframe.
#' @param min_sun numeric. minimum amount of sunlit pixels in the plots
#
#' @return list
#' @importFrom tidyr nest
#' @export
#'
spectralspecies_per_rast_sample <- function(input_rast, Kmeans_info,
                                            selected_bands = NULL,
                                            fd_metrics = NULL, input_mask = NULL,
                                            rast_sample = NULL,
                                            AttributeTable = NULL, min_sun = 0.25){

  FunctDiv <- SSValid <- NULL
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
    # Functional diversity
    if (!is.null(fd_metrics)){
      if (is.null(selected_bands))
        selected_bands <- seq_len(ncol(rast_sample_noID))
      # center reduce data
      inputdata_cr <- center_reduce(x = rast_sample_noID[selected_bands],
                                    m = Kmeans_info$MinVal,
                                    sig = Kmeans_info$Range)
      inputdata_cr$ID <- rast_sample$ID
      inputdata_cr <- inputdata_cr %>% split(.$ID)
      inputdata_cr <- lapply(inputdata_cr,
                             function(x) x[ , !(names(x) %in% "ID")])
      # in case only one dimension
      inputdata_cr <- lapply(inputdata_cr, data.frame)
      FunctDiv <- lapply(X = inputdata_cr,
                         FUN = get_functional_diversity,
                         fd_metrics = fd_metrics)
      FunctDiv <- data.frame('FRic' = unlist(lapply(FunctDiv, '[[',1)),
                             'FEve' = unlist(lapply(FunctDiv, '[[',2)),
                             'FDiv' = unlist(lapply(FunctDiv, '[[',3)),
                             'FDis' = unlist(lapply(FunctDiv, '[[',4)),
                             'FRaoq' = unlist(lapply(FunctDiv, '[[',5)))
    }
  }
  return(list('SSValid' = SSValid,
              'AttributeTable' = AttributeTable,
              'FunctDiv' = FunctDiv))
}
