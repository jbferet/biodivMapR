# ==============================================================================
# biodivMapR
# Lib_MapBetaDiversity.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
# ==============================================================================
# This Library computes bray curtis dissimilarity among spatial units based on
# spectral species distribution file and generates a RGB representation with NMDS
# ==============================================================================

#' maps beta diversity indicator based on spectral species distribution
#'
#' @param Input.Image.File character. Path and name of the image to be processed.
#' @param Output.Dir character. Output directory.
#' @param Spatial.Unit numeric. Dimensions of the spatial unit.
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...).
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param Nb.Units.Ordin numeric. Maximum number of spatial units to be processed in NMDS.
#' --> 1000 will be fast but may not capture important patterns if large area
#' --> 4000 will be slow but may show better ability to capture landscape patterns
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param scaling character. "PCO" or "NMDS".
#' @param FullRes boolean.
#' @param LowRes boolean.
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param pcelim numeric. Minimum contribution (in \%) required for a spectral species
#'
#' @export
map_beta_div <- function(Input.Image.File, Output.Dir, Spatial.Unit,
                               TypePCA = "SPCA", nbclusters = 50,
                               Nb.Units.Ordin = 2000, MinSun = 0.25,
                               pcelim = 0.02, scaling = "PCO", FullRes = TRUE,
                               LowRes = FALSE, nbCPU = 1, MaxRAM = 0.25) {
  Output.Dir2 <- Define.Output.Dir(Output.Dir, Input.Image.File, TypePCA)
  WS_Save <- paste(Output.Dir2, "PCA_Info.RData", sep = "")
  Output.Dir.SS <- Define.Output.SubDir(Output.Dir, Input.Image.File, TypePCA, "SpectralSpecies")
  Output.Dir.BETA <- Define.Output.SubDir(Output.Dir, Input.Image.File, TypePCA, "BETA")
  load(file = WS_Save)
  Beta <- Compute.BetaDiversity(Output.Dir.SS, MinSun, Nb.Units.Ordin, NbIter,
                                nbclusters, pcelim, scaling = scaling,
                                nbCPU = nbCPU, MaxRAM = MaxRAM)
  # Create images corresponding to Beta-diversity
  print("Write beta diversity maps")
  Index <- paste("BetaDiversity_BCdiss_", scaling, sep = "")
  Beta.Path <- paste(Output.Dir.BETA, Index, "_", Spatial.Unit, sep = "")
  Write.Image.Beta(Beta$BetaDiversity, Beta$HDR, Beta.Path, Spatial.Unit, FullRes = FullRes, LowRes = LowRes)
  return()
}

# maps beta diversity indicator based on spectral species distribution
# and saves in directories specifying the number of clusters
#
# @param Input.Image.File Path and name of the image to be processed
# @param Output.Dir output directory
# @param Spatial.Unit dimensions of the spatial unit
# @param TypePCA Type of PCA (PCA, SPCA, NLPCA...)
# @param nbclusters number of clusters defined in k-Means
# @param Nb.Units.Ordin maximum number of spatial units to be processed in NMDS
# @param MinSun minimum proportion of sunlit pixels required to consider plot
# @param pcelim minimum contribution (in %) required for a spectral species
# --> 1000 will be fast but may not capture important patterns if large area
# @param scaling
# @param FullRes
# @param LowRes
# @param nbCPU
# @param MaxRAM
# --> 4000 will be slow but may show better ability to capture landscape patterns
# @return
Map.Beta.Diversity.TestnbCluster <- function(Input.Image.File, Output.Dir, Spatial.Unit, TypePCA = "SPCA", nbclusters = 50, Nb.Units.Ordin = 2000, MinSun = 0.25, pcelim = 0.02, scaling = "PCO", FullRes = TRUE, LowRes = FALSE, nbCPU = FALSE, MaxRAM = FALSE) {
  Output.Dir2 <- Define.Output.Dir(Output.Dir, Input.Image.File, TypePCA)
  WS_Save <- paste(Output.Dir, "PCA_Info.RData", sep = "")
  Output.Dir.SS <- Define.Output.SubDir(Output.Dir, Input.Image.File, TypePCA, paste("SpectralSpecies_", nbclusters, sep = ""))
  Output.Dir.BETA <- Define.Output.SubDir(Output.Dir, Input.Image.File, TypePCA, paste("BETA_", nbclusters, sep = ""))
  load(file = WS_Save)
  Beta <- Compute.BetaDiversity(Output.Dir.SS, MinSun, Nb.Units.Ordin, NbIter, nbclusters, pcelim, scaling = scaling, nbCPU = nbCPU, MaxRAM = MaxRAM)
  # Create images corresponding to Beta-diversity
  print("Write beta diversity maps")
  Index <- paste("BetaDiversity_BCdiss_", scaling, sep = "")
  Beta.Path <- paste(Output.Dir.BETA, Index, "_", Spatial.Unit, sep = "")
  Write.Image.Beta(Beta$BetaDiversity, Beta$HDR, Beta.Path, Spatial.Unit, FullRes = FullRes, LowRes = LowRes)
  return()
}

# Gets sunlit pixels from SpectralSpecies_Distribution_Sunlit
#
# @param ImPathSunlit path for SpectralSpecies_Distribution_Sunlit
# @param MinSun minimum proportion of sunlit pixels required
#
# @return
Get.Sunlit.Pixels <- function(ImPathSunlit, MinSun) {

  # Filter out spatial units showing poor illumination
  Sunlit.HDR <- Get.HDR.Name(ImPathSunlit)
  HDR.Sunlit <- read.ENVI.header(Sunlit.HDR)
  nbpix <- as.double(HDR.Sunlit$lines) * as.double(HDR.Sunlit$samples)
  fid <- file(
    description = ImPathSunlit, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  Sunlit <- readBin(fid, numeric(), n = nbpix, size = 4)
  close(fid)
  Sunlit <- aperm(array(Sunlit, dim = c(HDR.Sunlit$samples, HDR.Sunlit$lines)))
  # define sunlit spatial units
  Select.Sunlit <- which(Sunlit > MinSun)
  # define where to extract each datapoint in the file
  coordi <- ((Select.Sunlit - 1) %% HDR.Sunlit$lines) + 1
  coordj <- floor((Select.Sunlit - 1) / HDR.Sunlit$lines) + 1
  # sort based on line and col (important for optimal scan of large files)
  coordTot <- cbind(coordi, coordj)
  # sort samples from top to bottom in order to optimize read/write of the image
  # indxTot saves the order of the data for reconstruction later
  indxTot <- order(coordTot[, 1], coordTot[, 2], na.last = FALSE)
  coordTotSort <- coordTot[indxTot, ]
  Select.Sunlit <- Select.Sunlit[indxTot]

  my_list <- list("Select.Sunlit" = Select.Sunlit, "coordTotSort" = coordTotSort)
  return(my_list)
}

# computes NMDS
#
# @param MatBCdist BC dissimilarity matrix
#
# @return BetaNMDS.sel
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
Compute.NMDS <- function(MatBCdist) {
  nbiterNMDS <- 4
  if (Sys.info()["sysname"] == "Windows") {
    nbCoresNMDS <- 2
  } else if (Sys.info()["sysname"] == "Linux") {
    nbCoresNMDS <- 4
  }
  # multiprocess of spectral species distribution and alpha diversity metrics
  plan(multiprocess, workers = nbCoresNMDS) ## Parallelize using four cores
  BetaNMDS <- future_lapply(MatBCdist, FUN = NMDS, mindim = 3, maxdim = 3, nits = 1, future.packages = c("ecodist"))
  plan(sequential)
  # find iteration with minimum stress
  Stress <- vector(length = nbiterNMDS)
  for (i in 1:nbiterNMDS) {
    Stress[i] <- BetaNMDS[[i]]$stress
  }
  print("Stress obtained for NMDS iterations:")
  print(Stress)
  print("Rule of thumb")
  print("stress < 0.05 provides an excellent represention in reduced dimensions")
  print("stress < 0.1 is great")
  print("stress < 0.2 is good")
  print("stress > 0.3 provides a poor representation")
  MinStress <- find(Stress == min(Stress))
  BetaNMDS.sel <- BetaNMDS[[MinStress]]$conf
  BetaNMDS.sel <- data.frame(BetaNMDS.sel[[1]])
  return(BetaNMDS.sel)
}

# identifies ordination coordinates based on nearest neighbors
#
# @param Beta.Ordination.sel ordination
# @param SSD.Path ath for spectral species distribution file
# @param Sample.Sel Samples selected during ordination
# @param coordTotSort coordinates of sunlit spatial units
# @param NbIter number of iterations
# @param nbclusters number of clusters
# @param pcelim number of CPUs available
# @param nbCPU number of CPUs available
#
# @return Ordination.est coordinates of each spatial unit in ordination space
#' @importFrom snow splitRows
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
Ordination.to.NN <- function(Beta.Ordination.sel, SSD.Path, Sample.Sel, coordTotSort, NbIter, nbclusters, pcelim, nbCPU = FALSE) {
  nb.Sunlit <- dim(coordTotSort)[1]
  # define number of samples to be sampled each time during paralle processing
  nb.samples.per.sub <- round(1e7 / dim(Sample.Sel)[1])
  # number of paralle processes to run
  nb.sub <- round(nb.Sunlit / nb.samples.per.sub)
  if (nb.sub == 0) nb.sub <- 1
  id.sub <- splitRows(as.matrix(seq(1, nb.Sunlit, by = 1), ncol = 1), ncl = nb.sub)
  # compute ordination coordinates from each subpart
  Nb.Units.Ordin <- dim(Sample.Sel)[1]
  plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
  Schedule.Per.Thread <- ceiling(nb.sub / nbCPU)
  OutPut <- future_lapply(id.sub,
    FUN = Ordination.Parallel, coordTotSort = coordTotSort, SSD.Path = SSD.Path,
    Sample.Sel = Sample.Sel, Beta.Ordination.sel = Beta.Ordination.sel, Nb.Units.Ordin = Nb.Units.Ordin,
    NbIter = NbIter, nbclusters = nbclusters, pcelim = pcelim, future.scheduling = Schedule.Per.Thread,
    future.packages = c("vegan", "dissUtils", "R.utils", "tools", "snow", "matlab")
  )
  plan(sequential)
  Ordination.est <- do.call("rbind", OutPut)
  gc()
  return(Ordination.est)
}

# applies results of ordination to full image based on nearest neighbors
#
# @param id.sub
# @param coordTotSort
# @param SSD.Path
# @param Sample.Sel
# @param Beta.Ordination.sel
# @param Nb.Units.Ordin
# @param NbIter
# @param nbclusters
# @param pcelim
#
# @return list of mean and SD of alpha diversity metrics
Ordination.Parallel <- function(id.sub, coordTotSort, SSD.Path, Sample.Sel, Beta.Ordination.sel, Nb.Units.Ordin, NbIter, nbclusters, pcelim) {

  # get Spectral species distribution
  coordPix <- coordTotSort[id.sub, ]
  SSD.NN <- Extract.Samples.From.Image(SSD.Path, coordPix)
  # compute the mean BC dissimilarity sequentially for each iteration
  MatBCtmp <- matrix(0, nrow = nrow(id.sub), ncol = Nb.Units.Ordin)
  SSDList <- list()
  for (i in 1:NbIter) {
    lb <- 1 + (i - 1) * nbclusters
    ub <- i * nbclusters
    SSDList[[1]] <- SSD.NN[, lb:ub]
    SSDList[[2]] <- Sample.Sel[, lb:ub]
    MatBCtmp0 <- Compute.BCdiss(SSDList, pcelim)
    MatBCtmp <- MatBCtmp + MatBCtmp0
  }
  MatBCtmp <- MatBCtmp / NbIter
  # get the knn closest neighbors from each kernel
  knn <- 3
  OutPut <- Compute.NN.From.Ordination(MatBCtmp, knn, Beta.Ordination.sel)
  return(OutPut)
}

# computes beta diversity
#
# @param Output.Dir directory where spectral species are stored
# @param MinSun minimum proportion of sunlit pixels required to consider plot
# @param Nb.Units.Ordin maximum number of spatial units to be processed in Ordination
# @param NbIter number of iterations
# @param nbclusters number of clusters defined in k-Means
# @param scaling
# @param nbCPU
# @param MaxRAM
# @param pcelim min proprtion for a spectral species in a spatial unit to be considerd
#
# @return
#' @importFrom labdsv pco
Compute.BetaDiversity <- function(Output.Dir, MinSun, Nb.Units.Ordin, NbIter, nbclusters, pcelim, scaling = "PCO", nbCPU = FALSE, MaxRAM = FALSE) {
  # Define path for images to be used
  SSD.Path <- paste(Output.Dir, "SpectralSpecies_Distribution", sep = "")
  ImPathSunlit <- paste(Output.Dir, "SpectralSpecies_Distribution_Sunlit", sep = "")
  # Get illuminated pixels based on  SpectralSpecies_Distribution_Sunlit
  Sunlit.Pixels <- Get.Sunlit.Pixels(ImPathSunlit, MinSun)
  Select.Sunlit <- Sunlit.Pixels$Select.Sunlit
  nb.Sunlit <- length(Select.Sunlit)
  # Define spatial units processed through ordination and those processed through
  # Nearest neighbor based on the first ordination batch
  print("Read Spectral Species distribution")
  RandPermKernels <- sample(seq(1, nb.Sunlit, by = 1))
  if (Nb.Units.Ordin <= nb.Sunlit) {
    Kernels.NN <- RandPermKernels[(Nb.Units.Ordin + 1):nb.Sunlit]
  } else {
    Nb.Units.Ordin <- nb.Sunlit
    Kernels.NN <- c()
  }
  # read spectral species distribution file
  SSD.All <- Extract.Samples.From.Image(SSD.Path, Sunlit.Pixels$coordTotSort)
  # define kernels used for Ordination
  Kernels.Ordination <- RandPermKernels[1:Nb.Units.Ordin]
  Sample.Sel <- SSD.All[Kernels.Ordination, ]
  rm(SSD.All)
  gc()

  # create a Bray curtis dissimilarity matrix for each iteration
  print("compute BC dissimilarity for selected kernels")
  # create a list in with each element is an iteration
  MatBC <- matrix(0, ncol = Nb.Units.Ordin, nrow = Nb.Units.Ordin)
  SSDList <- list()
  BC.from.SSD <- list()
  for (i in 1:NbIter) {
    lb <- 1 + (i - 1) * nbclusters
    ub <- i * nbclusters
    SSDList[[1]] <- Sample.Sel[, lb:ub]
    SSDList[[2]] <- Sample.Sel[, lb:ub]
    BC.from.SSD <- Compute.BCdiss(SSDList, pcelim)
    MatBC <- MatBC + BC.from.SSD
  }
  MatBC <- MatBC / NbIter

  # Perform Ordination based on BC dissimilarity matrix
  print("perform Ordination on the BC dissimilarity matrix averaged from all iterations")
  # parallel computing of Ordination can be run on 2 cores on Windows.
  # core management seems better on linux --> 4 cores possible
  MatBCdist <- as.dist(MatBC, diag = FALSE, upper = FALSE)
  if (scaling == "NMDS") {
    Beta.Ordination.sel <- Compute.NMDS(MatBCdist)
  } else if (scaling == "PCO") {
    BetaPCO <- pco(MatBCdist, k = 3)
    Beta.Ordination.sel <- BetaPCO$points
  }

  # Perform nearest neighbor on spatial units excluded from Ordination
  print("BC dissimilarity between samples selected for Ordination and remaining")
  coordTotSort <- Sunlit.Pixels$coordTotSort
  Ordination.est <- Ordination.to.NN(Beta.Ordination.sel, SSD.Path, Sample.Sel, coordTotSort, NbIter, nbclusters, pcelim, nbCPU = nbCPU)

  # Reconstuct spatialized beta diversity map from previous computation
  Sunlit.HDR <- Get.HDR.Name(ImPathSunlit)
  HDR.Sunlit <- read.ENVI.header(Sunlit.HDR)
  BetaDiversity <- as.matrix(Ordination.est, ncol = 3)
  BetaDiversityRGB <- array(NA, c(as.double(HDR.Sunlit$lines), as.double(HDR.Sunlit$samples), 3))
  BetaTmp <- matrix(NA, nrow = as.double(HDR.Sunlit$lines), ncol = as.double(HDR.Sunlit$samples))
  for (i in 1:3) {
    BetaTmp[Select.Sunlit] <- BetaDiversity[, i]
    BetaDiversityRGB[, , i] <- BetaTmp
  }
  list <- ls()
  rm(list = list[-which(list == "BetaDiversityRGB" | list == "Select.Sunlit" | list == "HDR.Sunlit")])
  gc()
  my_list <- list("BetaDiversity" = BetaDiversityRGB, "Select.Sunlit" = Select.Sunlit, "HDR" = HDR.Sunlit)
  return(my_list)
}

# compute the bray curtis dissimilarity matrix corresponding to a list of kernels
# (rows) defined by their spectral species (columns)
# SSDList is a list containing spectral species distribution for two sets of kernels ([[1]] and [[2]])
# pcelim is the threshold for minimum contributin of a spctral species to be kept
#
# @param SSDList list of 2 groups to compute BC dissimilarity from
# @param pcelim minimum proportion required for a species to be included (currently deactivated)
#
# @return MatBC matrix of bray curtis dissimilarity matrix
#' @importFrom dissUtils diss
Compute.BCdiss <- function(SSDList, pcelim) {
  # compute the proportion of each spectral species
  # Here, the proportion is computed with regards to the total number of sunlit pixels
  # One may want to determine if the results are similar when the proportion is computed
  # with regards to the total number of pixels (se*se)
  # however it would increase dissimilarity betwen kernels with different number of sunlit pixels
  SSD <- list()
  for (i in 1:length(SSDList)) {
    # get the total number of sunlit pixels in spatial unit
    SumSpecies <- rowSums(SSDList[[i]])
    # EDIT 02-Feb-2019: changed this
    # check if some kernels contain no data and eliminate if it is the case
    # elim          = which(SumSpecies==0)
    # if (length(elim)>0){
    #   SumSpecies  = SumSpecies[-elim]
    #   SSDList[[i]]= SSDList[[i]][-elim,]
    # }
    # SSD[[i]]      = apply(SSDList[[i]],2,function(x,c1) x/c1,"c1"=SumSpecies)
    # SSD[[i]][which(SSD[[i]]<pcelim)]  = 0
    elim <- which(SumSpecies == 0)
    if (length(elim) > 0) {
      SumSpecies[elim] <- 1
      SSDList[[i]][elim, ] <- 0
    }
    SSD[[i]] <- apply(SSDList[[i]], 2, function(x, c1) x / c1, "c1" = SumSpecies)
    SSD[[i]][which(SSD[[i]] < pcelim)] <- 0
  }
  # matrix of bray curtis dissimilarity (size = nb kernels x nb kernels)
  # Here use the package "dissUtils" to compute dissimilarity matrices sequentially
  MatBC <- diss(SSD[[1]], SSD[[2]], method = "braycurtis")
  # EDIT 06-Feb-2019: added this to fix problem when empty kernels occur, leading to NA BC value
  if (length(which(is.na(MatBC) == TRUE)) > 0) {
    MatBC[which(is.na(MatBC) == TRUE)] <- 1
  }
  return(MatBC)
}

# compute the nearest neighbors among kernels used in NMDS
#
# @param MatBC3 matrix of BC dissimilarity between the kernels excluded from Ordination (rows)
# @param knn number of neighbors
# @param BetaDiversity0
#
# @return estimated NMDS position based on nearest neighbors from NMDS
Compute.NN.From.Ordination <- function(MatBC3, knn, BetaDiversity0) {
  Ordin.est <- matrix(0, ncol = 3, nrow = nrow(MatBC3))
  for (i in 1:nrow(MatBC3)) {
    NNtmp <- sort(MatBC3[i, ], decreasing = FALSE, index.return = TRUE)
    Dist.Tot <- sum(NNtmp$x[1:knn])
    aa <- as.numeric(((Dist.Tot - NNtmp$x[1]) / (2 * Dist.Tot)) * BetaDiversity0[NNtmp$ix[1], ]
      + ((Dist.Tot - NNtmp$x[2]) / (2 * Dist.Tot)) * (BetaDiversity0[NNtmp$ix[2], ])
      + ((Dist.Tot - NNtmp$x[3]) / (2 * Dist.Tot)) * BetaDiversity0[NNtmp$ix[3], ])
    Ordin.est[i, 1:3] <- aa
  }
  return(Ordin.est)
}

# Writes image of beta diversity indicator (3 bands) resulting from BC + NMDS
#
# @param Image 2D matrix of image to be written
# @param HDR.SSD hdr template derived from SSD to modify
# @param ImagePath path of image file to be written
# @param Spatial.Unit spatial units use dto compute diversiy (in pixels)
# @param FullRes should full resolution image be written (original pixel size)
# @param LowRes should low resolution image be written (one value per spatial unit)
#
# @return
Write.Image.Beta <- function(Image, HDR.SSD, ImagePath, Spatial.Unit, FullRes = TRUE, LowRes = FALSE) {
  # Write image with resolution corresponding to Spatial.Unit
  HDR.Beta <- HDR.SSD
  HDR.Beta$bands <- 3
  HDR.Beta$`data type` <- 4
  PCs <- list()
  for (i in 1:3) {
    PCs <- c(PCs, paste("NMDS ", i))
  }
  PCs <- paste(PCs, collapse = ", ")
  HDR.Beta$`band names` <- PCs
  Image.Format <- ENVI.Type2Bytes(HDR.Beta)
  if (LowRes == TRUE) {
    headerFpath <- paste(ImagePath, ".hdr", sep = "")
    write.ENVI.header(HDR.Beta, headerFpath)
    ImgWrite <- aperm(Image, c(2, 3, 1))
    fidOUT <- file(
      description = ImagePath, open = "wb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    writeBin(c(ImgWrite), fidOUT, size = Image.Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
  }
  if (FullRes == TRUE) {
    # Write image with Full native resolution
    HDR.Full <- HDR.Beta
    HDR.Full$samples <- HDR.Beta$samples * Spatial.Unit
    HDR.Full$lines <- HDR.Beta$lines * Spatial.Unit
    HDR.Full <- Revert.Resolution.HDR(HDR.Full, Spatial.Unit)
    ImagePath.FullRes <- paste(ImagePath, "_Fullres", sep = "")
    headerFpath <- paste(ImagePath.FullRes, ".hdr", sep = "")
    write.ENVI.header(HDR.Full, headerFpath)
    Image.Format <- ENVI.Type2Bytes(HDR.Full)
    Image.FullRes <- array(NA, c(HDR.Full$lines, HDR.Full$samples, 3))
    for (b in 1:3) {
      for (i in 1:HDR.SSD$lines) {
        for (j in 1:HDR.SSD$samples) {
          Image.FullRes[((i - 1) * Spatial.Unit + 1):(i * Spatial.Unit), ((j - 1) * Spatial.Unit + 1):(j * Spatial.Unit), b] <- Image[i, j, b]
        }
      }
    }
    ImgWrite <- aperm(Image.FullRes, c(2, 3, 1))
    fidOUT <- file(
      description = ImagePath.FullRes, open = "wb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    writeBin(c(ImgWrite), fidOUT, size = Image.Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
    # zip resulting file
    ZipFile(ImagePath.FullRes)
  }
  return("")
}
