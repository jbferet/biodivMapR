# ==============================================================================
# biodivMapR
# Lib_MapBetaDiversity.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library computes bray curtis dissimilarity among spatial units based on
# spectral species distribution file and generates a RGB representation with NMDS
# ==============================================================================

#' maps beta diversity indicator based on spectral species distribution
#'
#' @param Input_Image_File character. Path and name of the image to be processed.
#' @param Output_Dir character. Output directory.
#' @param window_size numeric. Dimensions of the spatial unit.
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...).
#' @param nb_partitions numeric. Number of partitions (repetitions) to be computed then averaged.
#' @param nbclusters numeric. Number of clusters defined in k-Means.
#' @param Nb_Units_Ordin numeric. Maximum number of spatial units to be processed in NMDS.
#' --> 1000 will be fast but may not capture important patterns if large area
#' --> 4000 will be slow but may show better ability to capture landscape patterns
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param pcelim numeric. Minimum contribution (in \%) required for a spectral species
#' @param scaling character. "PCO" or "NMDS".
#' @param FullRes boolean.
#' @param LowRes boolean.
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#'
#' @return None
#' @export
map_beta_div <- function(Input_Image_File, Output_Dir, window_size,
                               TypePCA = "SPCA", nb_partitions = 20,nbclusters = 50,
                               Nb_Units_Ordin = 2000, MinSun = 0.25,
                               pcelim = 0.02, scaling = "PCO", FullRes = TRUE,
                               LowRes = FALSE, nbCPU = 1, MaxRAM = 0.25) {

  Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "SpectralSpecies")
  Output_Dir_BETA <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "BETA")
  Beta <- compute_beta_metrics(Output_Dir = Output_Dir_SS, MinSun = MinSun,
                               Nb_Units_Ordin = Nb_Units_Ordin, nb_partitions = nb_partitions,
                               nbclusters = nbclusters, pcelim = pcelim, scaling = scaling,
                               nbCPU = nbCPU, MaxRAM = MaxRAM)
  # Create images corresponding to Beta-diversity
  print("Write beta diversity maps")
  Index <- paste("BetaDiversity_BCdiss_", scaling, sep = "")
  Beta.Path <- paste(Output_Dir_BETA, Index, "_", window_size, sep = "")
  write_raster(Image = Beta$BetaDiversity, HDR = Beta$HDR, ImagePath = Beta.Path,
               window_size = window_size, FullRes = FullRes, LowRes = TRUE,
               SmoothImage = FALSE)
  return(invisible())
}

#' computes NMDS
#
#' @param MatBCdist BC dissimilarity matrix
#
#' @return BetaNMDS_sel
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom ecodist nmds
#' @importFrom utils find
compute_NMDS <- function(MatBCdist) {
  nbiterNMDS <- 4
  if (Sys.info()["sysname"] == "Windows") {
    nbCoresNMDS <- 2
  } else if (Sys.info()["sysname"] == "Linux") {
    nbCoresNMDS <- 4
  }
  # multiprocess of spectral species distribution and alpha diversity metrics
  plan(multiprocess, workers = nbCoresNMDS) ## Parallelize using four cores
  BetaNMDS <- future_lapply(MatBCdist, FUN = nmds, mindim = 3, maxdim = 3, nits = 1, future.packages = c("ecodist"))
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
  BetaNMDS_sel <- BetaNMDS[[MinStress]]$conf
  BetaNMDS_sel <- data.frame(BetaNMDS_sel[[1]])
  return(BetaNMDS_sel)
}

# identifies ordination coordinates based on nearest neighbors
#
# @param Beta_Ordination_sel ordination
# @param SSD_Path ath for spectral species distribution file
# @param Sample_Sel Samples selected during ordination
# @param coordTotSort coordinates of sunlit spatial units
# @param nb_partitions number of k-means then averaged
# @param nbclusters number of clusters
# @param pcelim number of CPUs available
# @param nbCPU number of CPUs available
#
# @return Ordination_est coordinates of each spatial unit in ordination space
#' @importFrom snow splitRows
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
ordination_to_NN <- function(Beta_Ordination_sel, SSD_Path, Sample_Sel, coordTotSort, nb_partitions, nbclusters, pcelim, nbCPU = FALSE) {
  nb_Sunlit <- dim(coordTotSort)[1]
  # define number of samples to be sampled each time during paralle processing
  nb_samples_per_sub <- round(1e7 / dim(Sample_Sel)[1])
  # number of paralle processes to run
  nb.sub <- round(nb_Sunlit / nb_samples_per_sub)
  if (nb.sub == 0) nb.sub <- 1
  id.sub <- splitRows(as.matrix(seq(1, nb_Sunlit, by = 1), ncol = 1), ncl = nb.sub)
  # compute ordination coordinates from each subpart
  Nb_Units_Ordin <- dim(Sample_Sel)[1]
  plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
  Schedule_Per_Thread <- ceiling(nb.sub / nbCPU)
  OutPut <- future_lapply(id.sub,
    FUN = ordination_parallel, coordTotSort = coordTotSort, SSD_Path = SSD_Path,
    Sample_Sel = Sample_Sel, Beta_Ordination_sel = Beta_Ordination_sel, Nb_Units_Ordin = Nb_Units_Ordin,
    nb_partitions = nb_partitions, nbclusters = nbclusters, pcelim = pcelim, future.scheduling = Schedule_Per_Thread,
    future.packages = c("vegan", "dissUtils", "R.utils", "tools", "snow", "matlab")
  )
  plan(sequential)
  Ordination_est <- do.call("rbind", OutPut)
  gc()
  return(Ordination_est)
}

# applies results of ordination to full image based on nearest neighbors
#
# @param id.sub
# @param coordTotSort
# @param SSD_Path
# @param Sample_Sel
# @param Beta_Ordination_sel
# @param Nb_Units_Ordin
# @param nb_partitions
# @param nbclusters
# @param pcelim
#
# @return list of mean and SD of alpha diversity metrics
ordination_parallel <- function(id.sub, coordTotSort, SSD_Path, Sample_Sel, Beta_Ordination_sel, Nb_Units_Ordin, nb_partitions, nbclusters, pcelim) {

  # get Spectral species distribution
  coordPix <- coordTotSort[id.sub, ]
  SSD_NN <- extract_samples_from_image(SSD_Path, coordPix)
  # compute the mean BC dissimilarity sequentially for each iteration
  MatBCtmp <- matrix(0, nrow = nrow(id.sub), ncol = Nb_Units_Ordin)
  SSDList <- list()
  for (i in 1:nb_partitions) {
    lb <- 1 + (i - 1) * nbclusters
    ub <- i * nbclusters
    SSDList[[1]] <- SSD_NN[, lb:ub]
    SSDList[[2]] <- Sample_Sel[, lb:ub]
    MatBCtmp0 <- compute_BCdiss(SSDList, pcelim)
    MatBCtmp <- MatBCtmp + MatBCtmp0
  }
  MatBCtmp <- MatBCtmp / nb_partitions
  # get the knn closest neighbors from each kernel
  knn <- 3
  OutPut <- compute_NN_from_ordination(MatBCtmp, knn, Beta_Ordination_sel)
  return(OutPut)
}

# computes beta diversity
#
# @param Output_Dir directory where spectral species are stored
# @param MinSun minimum proportion of sunlit pixels required to consider plot
# @param Nb_Units_Ordin maximum number of spatial units to be processed in Ordination
# @param nb_partitions number of iterations
# @param nbclusters number of clusters defined in k-Means
# @param scaling
# @param nbCPU
# @param MaxRAM
# @param pcelim min proprtion for a spectral species in a spatial unit to be considerd
#
# @return
#' @importFrom labdsv pco
#' @importFrom stats as.dist
compute_beta_metrics <- function(Output_Dir, MinSun, Nb_Units_Ordin, nb_partitions, nbclusters, pcelim, scaling = "PCO", nbCPU = FALSE, MaxRAM = FALSE) {
  # Define path for images to be used
  SSD_Path <- paste(Output_Dir, "SpectralSpecies_Distribution", sep = "")
  ImPathSunlit <- paste(Output_Dir, "SpectralSpecies_Distribution_Sunlit", sep = "")
  # Get illuminated pixels based on  SpectralSpecies_Distribution_Sunlit
  Sunlit_Pixels <- get_sunlit_pixels(ImPathSunlit, MinSun)
  Select_Sunlit <- Sunlit_Pixels$Select_Sunlit
  nb_Sunlit <- length(Select_Sunlit)
  # Define spatial units processed through ordination and those processed through
  # Nearest neighbor based on the first ordination batch
  print("Read Spectral Species distribution")
  RandPermKernels <- sample(seq(1, nb_Sunlit, by = 1))
  if (Nb_Units_Ordin <= nb_Sunlit) {
    Kernels_NN <- RandPermKernels[(Nb_Units_Ordin + 1):nb_Sunlit]
  } else {
    Nb_Units_Ordin <- nb_Sunlit
    Kernels_NN <- c()
  }
  # read spectral species distribution file
  SSD_All <- extract_samples_from_image(SSD_Path, Sunlit_Pixels$coordTotSort)
  # define kernels used for Ordination
  Kernels_Ordination <- RandPermKernels[1:Nb_Units_Ordin]
  Sample_Sel <- SSD_All[Kernels_Ordination, ]
  rm(SSD_All)
  gc()

  # create a Bray curtis dissimilarity matrix for each iteration
  print("compute BC dissimilarity for selected kernels")
  # create a list in with each element is an iteration
  MatBC <- matrix(0, ncol = Nb_Units_Ordin, nrow = Nb_Units_Ordin)
  SSDList <- list()
  BC.from.SSD <- list()
  for (i in 1:nb_partitions) {
    lb <- 1 + (i - 1) * nbclusters
    ub <- i * nbclusters
    SSDList[[1]] <- Sample_Sel[, lb:ub]
    SSDList[[2]] <- Sample_Sel[, lb:ub]
    BC.from.SSD <- compute_BCdiss(SSDList, pcelim)
    MatBC <- MatBC + BC.from.SSD
  }
  MatBC <- MatBC / nb_partitions

  # Perform Ordination based on BC dissimilarity matrix
  print("perform Ordination on the BC dissimilarity matrix averaged from all iterations")
  # parallel computing of Ordination can be run on 2 cores on Windows.
  # core management seems better on linux --> 4 cores possible
  MatBCdist <- as.dist(MatBC, diag = FALSE, upper = FALSE)
  if (scaling == "NMDS") {
    Beta_Ordination_sel <- compute_NMDS(MatBCdist)
    PCname <- 'NMDS'
  } else if (scaling == "PCO") {
    BetaPCO <- pco(MatBCdist, k = 3)
    Beta_Ordination_sel <- BetaPCO$points
    PCname <- 'PCoA'
  }

  # Perform nearest neighbor on spatial units excluded from Ordination
  print("BC dissimilarity between samples selected for Ordination and remaining")
  coordTotSort <- Sunlit_Pixels$coordTotSort
  Ordination_est <- ordination_to_NN(Beta_Ordination_sel, SSD_Path, Sample_Sel, coordTotSort, nb_partitions, nbclusters, pcelim, nbCPU = nbCPU)

  # Reconstuct spatialized beta diversity map from previous computation
  Sunlit_HDR <- get_HDR_name(ImPathSunlit)
  HDR_Sunlit <- read_ENVI_header(Sunlit_HDR)
  BetaDiversity <- as.matrix(Ordination_est, ncol = 3)
  BetaDiversityRGB <- array(NA, c(as.double(HDR_Sunlit$lines), as.double(HDR_Sunlit$samples), 3))
  BetaTmp <- matrix(NA, nrow = as.double(HDR_Sunlit$lines), ncol = as.double(HDR_Sunlit$samples))
  for (i in 1:3) {
    BetaTmp[Select_Sunlit] <- BetaDiversity[, i]
    BetaDiversityRGB[, , i] <- BetaTmp
  }
  list <- ls()

  # update HDR file
  HDR_Beta <- HDR_Sunlit
  HDR_Beta$bands <- 3
  HDR_Beta$`data type` <- 4
  PCs <- list()
  for (i in 1:3) {
    PCs <- c(PCs, paste(PCname,'#', i,sep = ''))
  }
  PCs <- paste(PCs, collapse = ", ")
  HDR_Beta$`band names` <- PCs

  rm(list = list[-which(list == "BetaDiversityRGB" | list == "Select_Sunlit" | list == "HDR_Beta")])
  gc()
  my_list <- list("BetaDiversity" = BetaDiversityRGB, "Select_Sunlit" = Select_Sunlit, "HDR" = HDR_Beta)
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
compute_BCdiss <- function(SSDList, pcelim) {
  # compute the proportion of each spectral species
  # Here, the proportion is computed with regards to the total number of sunlit pixels
  # One may want to determine if the results are similar when the proportion is computed
  # with regards to the total number of pixels (se*se)
  # however it would increase dissimilarity betwen kernels with different number of sunlit pixels
  SSD <- list()
  for (i in 1:length(SSDList)) {
    # get the total number of sunlit pixels in spatial unit
    SumSpecies <- rowSums(SSDList[[i]])
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
compute_NN_from_ordination <- function(MatBC3, knn, BetaDiversity0) {
  Ordin_est <- matrix(0, ncol = 3, nrow = nrow(MatBC3))
  for (i in 1:nrow(MatBC3)) {
    NNtmp <- sort(MatBC3[i, ], decreasing = FALSE, index.return = TRUE)
    Dist_Tot <- sum(NNtmp$x[1:knn])
    aa <- as.numeric(((Dist_Tot - NNtmp$x[1]) / (2 * Dist_Tot)) * BetaDiversity0[NNtmp$ix[1], ]
      + ((Dist_Tot - NNtmp$x[2]) / (2 * Dist_Tot)) * (BetaDiversity0[NNtmp$ix[2], ])
      + ((Dist_Tot - NNtmp$x[3]) / (2 * Dist_Tot)) * BetaDiversity0[NNtmp$ix[3], ])
    Ordin_est[i, 1:3] <- aa
  }
  return(Ordin_est)
}

# Gets sunlit pixels from SpectralSpecies_Distribution_Sunlit
#
# @param ImPathSunlit path for SpectralSpecies_Distribution_Sunlit
# @param MinSun minimum proportion of sunlit pixels required
#
# @return
get_sunlit_pixels <- function(ImPathSunlit, MinSun) {

  # Filter out spatial units showing poor illumination
  Sunlit_HDR <- get_HDR_name(ImPathSunlit)
  HDR_Sunlit <- read_ENVI_header(Sunlit_HDR)
  nbpix <- as.double(HDR_Sunlit$lines) * as.double(HDR_Sunlit$samples)
  fid <- file(
    description = ImPathSunlit, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  Sunlit <- readBin(fid, numeric(), n = nbpix, size = 4)
  close(fid)
  Sunlit <- aperm(array(Sunlit, dim = c(HDR_Sunlit$samples, HDR_Sunlit$lines)))
  # define sunlit spatial units
  Select_Sunlit <- which(Sunlit > MinSun)
  # define where to extract each datapoint in the file
  coordi <- ((Select_Sunlit - 1) %% HDR_Sunlit$lines) + 1
  coordj <- floor((Select_Sunlit - 1) / HDR_Sunlit$lines) + 1
  # sort based on line and col (important for optimal scan of large files)
  coordTot <- cbind(coordi, coordj)
  # sort samples from top to bottom in order to optimize read/write of the image
  # indxTot saves the order of the data for reconstruction later
  indxTot <- order(coordTot[, 1], coordTot[, 2], na.last = FALSE)
  coordTotSort <- coordTot[indxTot, ]
  Select_Sunlit <- Select_Sunlit[indxTot]

  my_list <- list("Select_Sunlit" = Select_Sunlit, "coordTotSort" = coordTotSort)
  return(my_list)
}
