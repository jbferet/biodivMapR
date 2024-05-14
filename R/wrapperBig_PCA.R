#' this function aims at applying PCA on a raster or list of rasters
#' in combination with the function apply_bigRaster
#'
#' @param input_data list. Image data chunk and corresponding mask if available
#' @param input_args list. PCA model
#'
#' @return list. PCA transformed spectral information
#' @export
#'
wrapperBig_PCA<- function(input_data, input_args){

  CR <- input_args$CR
  Spectral <- input_args$Spectral
  BandsNoVar <- Spectral$BandsNoVar
  Bands2Keep <- Spectral$Bands2Keep
  PCA_model <- input_args$PCA_model
  Nb_PCs <- input_args$Nb_PCs
  nbPix <- dim(input_data[[1]])[1]
  output_data <- list('PCA' = matrix(NA, ncol = Nb_PCs, nrow = nbPix))
  # raster list from which PCA is computed
  nameBands <- names(input_data)
  if (!is.null(input_data$mask)) nameBands <- nameBands[-which(nameBands=='mask')]
  # discard masked pixels if necessary
  SelectPixels <- seq_len(nbPix)
  if (!is.null(input_data$mask)) SelectPixels <- which(input_data$mask>0)
  input_data$mask <- NULL
  if (length(SelectPixels)>0){
    input_data[[1]] <- input_data[[1]][SelectPixels,]
    # select spectral bands of interest
    if (CR){
      input_data[[1]] <- apply_continuum_removal(Spectral_Data = input_data[[1]],
                                                 Spectral = Spectral)
    } else {
      input_data[[1]] <- input_data[[1]][[Bands2Keep]]
    }
    if (length(BandsNoVar)>0) input_data[[1]] <- input_data[[1]][[-BandsNoVar]]
    PC_data <- scale(x = input_data[[1]],
                     center = PCA_model$center,
                     scale = PCA_model$scale) %*% PCA_model$rotation[, 1:Nb_PCs]
    output_data$PCA[SelectPixels, ] <- PC_data
  }
  rm(list=setdiff(ls(), "output_data"));gc()
  return(output_data)
}
