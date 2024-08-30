#' Computes continuum removal for matrix shaped data: more efficient than
#' processing individual spectra
#' the convex hull is based on the computation of the derivative between R at a
#' given spectral band and R at the following bands
#'
#' @param Minit numeric. initial data matrix (nb samples x nb bands)
#' @param Spectral_Bands numeric. central wavelength for the spectral bands
#' @param p list. progressor object for progress bar
#
#' @return samples from image and updated number of pixels to sampel if necessary
#' @export

continuumRemoval <- function(Minit, Spectral_Bands, p = NULL) {

  # Filter and prepare data prior to continuum removal
  CR_data <- filter_prior_CR(Minit, Spectral_Bands)
  Minit <- CR_data$Minit
  nbBands <- dim(Minit)[2]
  CR_data$Minit <- c()
  Spectral_Bands <- CR_data$Spectral_Bands
  nbSamples <- CR_data$nbSamples
  nbSamplesUpDate <- length(CR_data$SamplesToKeep)
  # if samples to be considered
  if (nbSamples > 0) {
    # initialization:
    # spectral band corresponding to each element of the data matrix
    Lambda <- repmat(matrix(Spectral_Bands, nrow = 1), nbSamplesUpDate, 1)
    # prepare matrices used to check evolution of the CR process
    # - elements still not processed through continuum removal: initialization to 1
    Still.Need.CR <- matrix(1, nrow = nbSamplesUpDate, ncol = nbBands)
    # - value of the convex hull: initially set to 0
    Convex_Hull <- matrix(0, nrow = nbSamplesUpDate, ncol = nbBands)
    # - reflectance value for latest interception with convex hull:
    # initialization to value of the first reflectance measurement
    Intercept_Hull <- repmat(matrix(Minit[, 1], ncol = 1), 1, nbBands)
    # - spectral band of latest interception
    Latest.Intercept <- repmat(X = matrix(Spectral_Bands[1], ncol = 1),
                               m = nbSamplesUpDate, n = nbBands)
    # number of spectral bands found as intercept
    nb.Intercept <- 0
    # continues until arbitrary stopping criterion:
    # stops when reach last spectral band (all values before last = 0)
    # while (max(Still.Need.CR[, seq_len(nbBands - 2)]) == 1 & (nb.Intercept <= (nbBands / 2))) {
    while (max(Still.Need.CR[, seq_len((nbBands - 2))]) == 1) {
      nb.Intercept <- nb.Intercept + 1
      # identify samples still needing continuum removal
      Sel <- which(Still.Need.CR[,(nbBands-2)]==1)
      # update variables to process samples needing CR only
      nbSamplesUpDate_tmp <- length(Sel)
      Lambda_tmp <- Lambda[Sel,]
      Minit_tmp <- Minit[Sel,]
      Latest.Intercept_tmp <- Latest.Intercept[Sel,]
      Still.Need.CR_tmp <- Still.Need.CR[Sel,]
      Convex_Hull_tmp <- Convex_Hull[Sel,]
      Intercept_Hull_tmp <- Intercept_Hull[Sel,]
      # Mstep give the position of the values to be updated
      Update_Data <- matrix(1, nrow = nbSamplesUpDate_tmp, ncol = nbBands)
      Update_Data[, nbBands] <- 0
      # initial step: first column set to 0; following steps: all bands below
      # max of the convex hull are set to 0
      Update_Data[which((Lambda_tmp - Latest.Intercept_tmp) < 0)] <- 0
      # compute slope for each coordinate
      Slope <- as.matrix((Minit_tmp - Intercept_Hull_tmp) / (Lambda_tmp - Latest.Intercept_tmp) * Still.Need.CR_tmp)
      # set current spectral band and previous bands to -9999
      if (!length(which(Still.Need.CR_tmp == 0)) == 0) {
        Slope[which(Still.Need.CR_tmp == 0)] <- -9999
      }
      if (!length(which(is.na(Slope))) == 0) {
        Slope[which(is.na(Slope))] <- -9999
      }
      # get max index for each row and convert into linear index
      Index.Max.Slope <- RowToLinear(max.col(Slope, ties.method = "last"),
                                     nbSamplesUpDate_tmp, nbBands)
      # !!!! OPTIM: replace repmat with column operation
      # update coordinates of latest intercept
      Latest.Intercept_tmp <- repmat(matrix(Lambda_tmp[Index.Max.Slope], ncol = 1), 1, nbBands)
      # update latest intercept
      Intercept_Hull_tmp <- repmat(matrix(as.matrix(Minit_tmp)[Index.Max.Slope], ncol = 1), 1, nbBands)
      # values corresponding to the domain between the two continuum maxima
      Update_Data[which((Lambda_tmp - Latest.Intercept_tmp) >= 0 | Latest.Intercept_tmp == Spectral_Bands[nbBands])] <- 0
      # values to eliminate for the next analysis: all spectral bands before latest intercept
      Still.Need.CR_tmp[which((Lambda_tmp - Latest.Intercept_tmp) < 0)] <- 0
      # the max slope is known, as well as the coordinates of the beginning and ending
      # a matrix now has to be built
      Convex_Hull_tmp <- Convex_Hull_tmp + Update_Data * (Intercept_Hull_tmp + sweep((Lambda_tmp - Latest.Intercept_tmp), 1, Slope[Index.Max.Slope], "*"))
      # update variables
      Convex_Hull[Sel,] <- Convex_Hull_tmp
      Still.Need.CR[Sel,] <- Still.Need.CR_tmp
      Lambda[Sel,] <- Lambda_tmp
      Latest.Intercept[Sel,] <- Latest.Intercept_tmp
      Intercept_Hull[Sel,] <- Intercept_Hull_tmp
    }
    CR_Results0 <- Minit[, 2:(nbBands - 2)] / Convex_Hull[, 2:(nbBands - 2)]
    CR_Results <- matrix(0, ncol = (nbBands - 3), nrow = nbSamples)
    CR_Results[CR_data$SamplesToKeep, ] <- as.matrix(CR_Results0)
  } else {
    CR_Results <- matrix(0, ncol = (nbBands - 3), nrow = nbSamples)
  }
  CR_Results <- data.frame(CR_Results)
  if (!is.null(p)) p()
  list <- ls()
  rm(list = list[-which(list == "CR_Results")])
  rm(list)
  gc()
  return(CR_Results)
}
