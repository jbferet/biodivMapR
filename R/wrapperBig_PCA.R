#' this function aims at applying PCA on a raster or list of rasters
#' in combination with the function apply_bigRaster
#'
#' @param input_data list. Image data chunk and corresponding mask if available
#' @param input_args list. PCA model and associated parameters required by bigRaster
#'
#' @return list. PCA transformed spectral information
#' @export
#'
wrapperBig_PCA<- function(input_data, input_args){

  spectral <- input_args$Spectral
  BandsNoVar <- spectral$BandsNoVar
  Bands2Keep <- spectral$Bands2Keep
  PCA_model <- input_args$PCA_model
  Nb_PCs <- input_args$Nb_PCs
  nb_pix <- dim(input_data[[1]])[1]
  output_data <- list('PCA' = matrix(NA, ncol = Nb_PCs, nrow = nb_pix))
  # raster list from which PCA is computed
  nameBands <- names(input_data)
  if (!is.null(input_data$mask))
    nameBands <- nameBands[-which(nameBands=='mask')]
  # discard masked pixels if necessary
  select_pixels <- seq_len(nb_pix)
  if (!is.null(input_data$mask))
    select_pixels <- which(input_data$mask>0)
  input_data$mask <- NULL
  if (length(select_pixels)>0){
    input_data[[1]] <- input_data[[1]][select_pixels,]
    # select spectral bands of interest
    if (input_args$CR){
      input_data[[1]] <- apply_continuum_removal(spectral_data = input_data[[1]],
                                                 spectral = spectral)
    } else {
      input_data[[1]] <- input_data[[1]][Bands2Keep]
    }
    if (length(BandsNoVar)>0)
      input_data[[1]] <- input_data[[1]][-BandsNoVar]
    PC_data <- scale(x = input_data[[1]],
                     center = PCA_model$center,
                     scale = PCA_model$scale) %*%
      PCA_model$rotation[, seq_len(Nb_PCs)]
    output_data$PCA[select_pixels, ] <- PC_data
  }
  rm(list=setdiff(ls(), "output_data"));gc()
  return(output_data)
}
