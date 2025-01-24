#' Performs radiometric filtering based on three criteria: NDVI, NIR reflectance, Blue reflectance
#'
#' @param blk character. Path of the image to be processed
#' @param r_in character. Path of the mask corresponding to the image
#' @param Thresholds character. Path for output directory
#' @param Spectral_Bands character. Type of PCA: choose either "PCA" or "SPCA"
#'
#' @return mask = updated mask values
#' @importFrom terra readValues
#' @export

radiometricfilter_chunk <- function(blk, r_in, Thresholds, Spectral_Bands){
  # 1- read input files
  input_data <- list()
  nameVars <- c()
  ll <- names(r_in)
  for (fid in ll) input_data[[fid]] <- terra::readValues(r_in[[fid]],
                                                         row = blk$row,
                                                         nrows = blk$nrows,
                                                         dataframe = TRUE)[[1]]
  names(input_data) <- ll
  mask <- 0*input_data[[1]] + 1
  for (band in names(Thresholds)){
    if (band == 'Blue '){
      elim <- which(input_data[[band]]>Thresholds[[band]])
      if (length(elim)>0) mask[elim] <- 0
    }
    if (band == 'NIR'){
      elim <- which(input_data[[band]]<Thresholds[[band]])
      if (length(elim)>0) mask[elim] <- 0
    }
    if (band == 'NDVI'){
      if (!is.null(input_data$NIR) & !is.null(input_data$Red)){
        NDVI <- (input_data$NIR-input_data$Red)/(input_data$NIR+input_data$Red)
        elim <- which(NDVI<Thresholds[[band]])
        if (length(elim)>0) mask[elim] <- 0
      }
    }
  }
  if (!is.null(input_data$mask)){
    elim <- which(input_data$mask==0)
    if (length(elim)>0) mask[elim] <- 0
  }
  return(mask)
}
