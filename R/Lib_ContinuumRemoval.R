# ===============================================================================
# biodivMapR
# Lib_ContinuumRemoval.R
# ===============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
# ===============================================================================
# This Library is dedicated to the computation of the continuum removal
# ===============================================================================

# prepares data to run multithreaded continuum removal
#
# @param Spectral.Data initial data matrix (nb samples x nb bands)
# @param nbCPU
# @param Spectral information about spectral bands
#
# @return samples from image and updated number of pixels to sampel if necessary
#' @importFrom snow splitRows
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
Apply.Continuum.Removal <- function(Spectral.Data, Spectral, nbCPU = 1) {
  if (!length(Spectral$WaterVapor) == 0) {
    Spectral.Data <- Spectral.Data[, -Spectral$WaterVapor]
  }

  # split data to perform continuum removal on into reasonable amount of data
  nb.Values <- dim(Spectral.Data)[1] * dim(Spectral.Data)[2]
  if (nb.Values > 0) {
    # corresponds to ~ 40 Mb data, but CR tends to requires ~ 10 times memory
    # avoids memory crash
    Max.nb.Values <- 2e6
    nb.CR <- ceiling(nb.Values / Max.nb.Values)
    Spectral.Data <- splitRows(Spectral.Data, nb.CR)
    # perform multithread continuum removal
    plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
    Schedule.Per.Thread <- ceiling(nb.CR / nbCPU)
    Spectral.Data.tmp <- future_lapply(Spectral.Data, FUN = ContinuumRemoval, Spectral.Bands = Spectral$Wavelength, future.scheduling = Schedule.Per.Thread)
    plan(sequential)
    Spectral.Data <- do.call("rbind", Spectral.Data.tmp)
    rm(Spectral.Data.tmp)
  } else {
    # edit 31-jan-2018
    # resize to delete first and last band as in continuum removal
    Spectral.Data <- Spectral.Data[, -c(1, 2)]
  }
  gc()
  return(Spectral.Data)
}

# Computes continuum removal for matrix shaped data: more efficient than
# processing individual spectra
# the convex hull is based on the computation of the derivative between R at a
# given spectral band and R at the following bands
#
# @param Minit initial data matrix (nb samples x nb bands)
# @param Spectral.Bands information about spectral bands
#
# @return samples from image and updated number of pixels to sampel if necessary
ContinuumRemoval <- function(Minit, Spectral.Bands) {

  # Filter and prepare data prior to continuum removal
  CR.data <- Filter.Prior.CR(Minit, Spectral.Bands)
  Minit <- CR.data$Minit
  nb.Bands <- dim(Minit)[2]
  CR.data$Minit <- c()
  Spectral.Bands <- CR.data$Spectral.Bands
  nb.Samples <- CR.data$nb.Samples
  # if samples to be considered
  if (nb.Samples > 0) {
    # initialization:
    # spectral band corresponding to each element of the data matrix
    Lambda <- repmat(matrix(Spectral.Bands, nrow = 1), nb.Samples, 1)
    # prepare matrices used to check evolution of the CR process
    # - elements still not processed through continuum removal: initialization to 1
    Still.Need.CR <- matrix(1, nrow = nb.Samples, ncol = nb.Bands)
    # - value of the convex hull: initially set to 0
    Convex.Hull <- matrix(0, nrow = nb.Samples, ncol = nb.Bands)
    # - reflectance value for latest interception with convex hull:
    # initialization to value of the first reflectance measurement
    Intercept.Hull <- repmat(matrix(Minit[, 1], ncol = 1), 1, nb.Bands)
    # - spectral band of latest interception
    Latest.Intercept <- repmat(matrix(Spectral.Bands[1], ncol = 1), nb.Samples, nb.Bands)
    # number of spectral bands found as intercept
    nb.Intercept <- 0
    # continues until arbitrary stopping criterion:
    # stops when reach last spectral band (all values before last = 0)
    while (max(Still.Need.CR[, 1:(nb.Bands - 2)]) == 1 & (nb.Intercept <= (nb.Bands / 2))) {
      nb.Intercept <- nb.Intercept + 1
      # Mstep give the position of the values to be updated
      Update.Data <- matrix(1, nrow = nb.Samples, ncol = nb.Bands)
      Update.Data[, nb.Bands] <- 0
      # initial step: first column set to 0; following steps: all bands below
      # max of the convex hull are set to 0
      Update.Data[which((Lambda - Latest.Intercept) < 0)] <- 0
      # compute slope for each coordinate
      Slope <- (Minit - Intercept.Hull) / (Lambda - Latest.Intercept) * Still.Need.CR
      # set current spectral band and previous bands to -9999
      if (!length(which(Still.Need.CR == 0)) == 0) {
        Slope[which(Still.Need.CR == 0)] <- -9999
      }
      if (!length(which(is.na(Slope))) == 0) {
        Slope[which(is.na(Slope))] <- -9999
      }
      # get max index for each row and convert into linear index
      Index.Max.Slope <- RowToLinear(max.col(Slope, ties.method = "last"), nb.Samples, nb.Bands)
      # !!!! OPTIM: replace repmat with column operation
      # update coordinates of latest intercept
      Latest.Intercept <- repmat(matrix(Lambda[Index.Max.Slope], ncol = 1), 1, nb.Bands)
      # update latest intercept
      Intercept.Hull <- repmat(matrix(Minit[Index.Max.Slope], ncol = 1), 1, nb.Bands)
      # values corresponding to the domain between the two continuum maxima
      Update.Data[which((Lambda - Latest.Intercept) >= 0 | Latest.Intercept == Spectral.Bands[nb.Bands])] <- 0
      # values to eliminate for the next analysis: all spectral bands before latest intercept
      Still.Need.CR[which((Lambda - Latest.Intercept) < 0)] <- 0
      # the max slope is known, as well as the coordinates of the beginning and ending
      # a matrix now has to be built
      Convex.Hull <- Convex.Hull + Update.Data * (Intercept.Hull + sweep((Lambda - Latest.Intercept), 1, Slope[Index.Max.Slope], "*"))
    }
    CR_Results0 <- Minit[, 2:(nb.Bands - 2)] / Convex.Hull[, 2:(nb.Bands - 2)]
    CR_Results <- matrix(0, ncol = (nb.Bands - 3), nrow = nb.Samples)
    CR_Results[CR.data$Samples.To.Keep, ] <- CR_Results0
  } else {
    CR_Results <- matrix(0, ncol = (nb.Bands - 3), nrow = nb.Samples)
  }
  list <- ls()
  rm(list = list[-which(list == "CR_Results")])
  rm(list)
  gc()
  return(CR_Results)
}

# Filter data prior to continuum removal:
# - values are expected to be real reflectance values between 0 and 10000
# - negative values may occur, so a +100 value is applied to avoid negative
# - possibly remaining negative values are set to 0
# - constant spectra are eliminated
#
# @param Spectral.Bands
# @param Minit initial data matrix, n rows = n samples, p cols = p spectral bands
#
# @return updated Minit
#' @importFrom matrixStats rowSds
Filter.Prior.CR <- function(Minit, Spectral.Bands) {

  # number of samples to be processed
  nb.Samples <- nrow(Minit)
  # make sure there is no negative values
  Minit <- Minit + 100.0
  Minit[which(Minit < 0)] <- 0
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
  # add negative values to the last column and update spectral bands
  Minit <- cbind(Minit, matrix(-9999, ncol = 1, nrow = nb.Samples))
  nb.Bands <- ncol(Minit)
  Spectral.Bands[nb.Bands] <- Spectral.Bands[nb.Bands - 1] + 10
  my_list <- list("Minit" = Minit, "Spectral.Bands" = Spectral.Bands, "nb.Samples" = nb.Samples, "Samples.To.Keep" = keep)
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
