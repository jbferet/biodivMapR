# ===============================================================================
# biodivMapR
# Lib_ContinuumRemoval.R
# ===============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/06 Jean-Baptiste FERET
# ===============================================================================
# This Library is dedicated to the computation of the continuum removal
# ===============================================================================

#' prepares data to run multithreaded continuum removal
#'
#' @param Spectral_Data numeric. initial data matrix (nb samples x nb bands)
#' @param Spectral list. information about spectral bands
#' @param nbCPU numeric. number of CPUs to be used in parallel
#
#' @return samples from image and updated number of pixels to sample if necessary
#' @importFrom snow splitRows
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor handlers with_progress
#' @export

apply_continuum_removal <- function(Spectral_Data, Spectral, nbCPU = 1) {
  if (!length(Spectral$WaterVapor) == 0) {
    Spectral_Data <- Spectral_Data[, -Spectral$WaterVapor]
  }

  # split data to perform continuum removal on into reasonable amount of data
  nb.Values <- dim(Spectral_Data)[1] * dim(Spectral_Data)[2]
  if (nb.Values > 0) {
    # corresponds to ~ 40 Mb data, but CR tends to requires ~ 10 times memory
    # avoids memory crash
    Max.nb.Values <- 2e6
    nb_CR <- ceiling(nb.Values / Max.nb.Values)
    Spectral_Data <- snow::splitRows(Spectral_Data, nb_CR)
    # perform multithread continuum removal
    if (nbCPU>1){
      future::plan(multisession, workers = nbCPU)
      handlers(global = TRUE)
      handlers("cli")
      with_progress({
        p <- progressr::progressor(steps = nb_CR)
        Spectral_Data_tmp <- future.apply::future_lapply(Spectral_Data,
                                                         FUN = continuumRemoval,
                                                         Spectral_Bands = Spectral$Wavelength,
                                                         p = p)
      })
      future::plan(sequential)
    } else {
      Spectral_Data_tmp <- lapply(Spectral_Data, FUN = continuumRemoval,
                                  Spectral_Bands = Spectral$Wavelength)
    }
    Spectral_Data <- do.call("rbind", Spectral_Data_tmp)
    rm(Spectral_Data_tmp)
  } else {
    # edit 31-jan-2018
    # resize to delete first and last band as in continuum removal
    Spectral_Data <- Spectral_Data[, -c(1, 2)]
  }
  gc()
  return(Spectral_Data)
}

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
    Latest.Intercept <- repmat(matrix(Spectral_Bands[1], ncol = 1), nbSamplesUpDate, nbBands)
    # number of spectral bands found as intercept
    nb.Intercept <- 0
    # continues until arbitrary stopping criterion:
    # stops when reach last spectral band (all values before last = 0)
    # while (max(Still.Need.CR[, 1:(nbBands - 2)]) == 1 & (nb.Intercept <= (nbBands / 2))) {
    while (max(Still.Need.CR[, 1:(nbBands - 2)]) == 1) {
      nb.Intercept <- nb.Intercept + 1
      # identify samples still needing continuum removal
      Sel <- which(Still.Need.CR[,(nbBands-2)]==1)
      # update variables to process samples needing CR only
      nbSamplesUpDate_tmp <- length(Sel)
      Lambda_tmp <- matrix(Lambda[Sel,],nrow = nbSamplesUpDate_tmp)
      Minit_tmp <- matrix(Minit[Sel,],nrow = nbSamplesUpDate_tmp)
      Latest.Intercept_tmp <- matrix(Latest.Intercept[Sel,],nrow = nbSamplesUpDate_tmp)
      Still.Need.CR_tmp <-matrix(Still.Need.CR[Sel,],nrow = nbSamplesUpDate_tmp)
      Convex_Hull_tmp <-matrix(Convex_Hull[Sel,],nrow = nbSamplesUpDate_tmp)
      Intercept_Hull_tmp <-matrix(Intercept_Hull[Sel,],nrow = nbSamplesUpDate_tmp)
      # Mstep give the position of the values to be updated
      Update_Data <- matrix(1, nrow = nbSamplesUpDate_tmp, ncol = nbBands)
      Update_Data[, nbBands] <- 0
      # initial step: first column set to 0; following steps: all bands below
      # max of the convex hull are set to 0
      Update_Data[which((Lambda_tmp - Latest.Intercept_tmp) < 0)] <- 0
      # compute slope for each coordinate
      Slope <- (Minit_tmp - Intercept_Hull_tmp) / (Lambda_tmp - Latest.Intercept_tmp) * Still.Need.CR_tmp
      # set current spectral band and previous bands to -9999
      if (!length(which(Still.Need.CR_tmp == 0)) == 0) {
        Slope[which(Still.Need.CR_tmp == 0)] <- -9999
      }
      if (!length(which(is.na(Slope))) == 0) {
        Slope[which(is.na(Slope))] <- -9999
      }
      # get max index for each row and convert into linear index
      Index.Max.Slope <- RowToLinear(max.col(Slope, ties.method = "last"), nbSamplesUpDate_tmp, nbBands)
      # !!!! OPTIM: replace repmat with column operation
      # update coordinates of latest intercept
      Latest.Intercept_tmp <- repmat(matrix(Lambda_tmp[Index.Max.Slope], ncol = 1), 1, nbBands)
      # update latest intercept
      Intercept_Hull_tmp <- repmat(matrix(Minit_tmp[Index.Max.Slope], ncol = 1), 1, nbBands)
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
    CR_Results[CR_data$SamplesToKeep, ] <- CR_Results0
  } else {
    CR_Results <- matrix(0, ncol = (nbBands - 3), nrow = nbSamples)
  }
  if (!is.null(p)){p()}
  list <- ls()
  rm(list = list[-which(list == "CR_Results")])
  rm(list)
  gc()
  return(CR_Results)
}

#' Filter data prior to continuum removal:
#' - values are expected to be real reflectance values between 0 and 10000
#' - negative values may occur, so a +100 value is applied to avoid negative
#' - possibly remaining negative values are set to 0
#' - constant spectra are eliminated
#'
#' @param Minit initial data matrix, n rows = n samples, p cols = p spectral bands
#' @param Spectral_Bands numeric. central wavelength for the spectral bands
#
#' @return list. updated Minit
#' @importFrom matrixStats rowSds
#' @export

filter_prior_CR <- function(Minit, Spectral_Bands) {

  # number of samples to be processed
  nbSamples <- nrow(Minit)
  # make sure there is no negative values
  # Minit[Minit<0] <- Minit + 100.0
  Minit[Minit < 0] <- 0
  # eliminate invariant spectra
  SD <- rowSds(Minit)
  elim <- which(SD == 0 | is.na(SD))
  keep <- which(!SD == 0 & !is.na(SD))
  nbi0 <- nrow(Minit)
  nbj <- ncol(Minit)
  if (nbj == 1 & nbi0 > 1) {
    Minit <- matrix(Minit, nrow = 1)
    nbi0 <- 1
    nbj <- ncol(Minit)
  }
  Minit <- Minit[keep, ]
  if (length(keep) == 1) {
    Minit <- matrix(Minit, nrow = 1)
  }
  nbSamplesUpDate <- nrow(Minit)
  # add negative values to the last column and update spectral bands
  Minit <- cbind(Minit, matrix(-9999, ncol = 1, nrow = nbSamplesUpDate))
  nbBands <- ncol(Minit)
  Spectral_Bands[nbBands] <- Spectral_Bands[nbBands - 1] + 10
  my_list <- list("Minit" = Minit, "Spectral_Bands" = Spectral_Bands, "nbSamples" = nbSamples, "SamplesToKeep" = keep)
  return(my_list)
}

#
# @param MM
# @param nbi
# @param nbj
#
# @return
RowToLinear <- function(MM, nbi, nbj) {
  adj <- seq(1:nbi)
  MaxCont <- ((MM - 1) * (nbi)) + adj
  return(MaxCont)
}

# R equivalent of repmat (matlab)
#
# @param X initial matrix
# @param m nb of replications in row dimension
# @param n nb of replications in column dimension
#
# @return
repmat <- function(X, m, n) {
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = T)
}
