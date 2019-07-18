# ==============================================================================
# biodivMapR
# Lib_MapSpectralSpecies.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2019/06 Jean-Baptiste FERET
# ==============================================================================
# This Library applies clustering on a selection of components stored in a PCA
# file previously created ("Perform_PCA.R") and produces corresponding spectral
# species
# ==============================================================================

#' maps spectral species based on PCA file computed previously
#'
#' @param Input.Image.File Path and name of the image to be processed
#' @param Output.Dir output directory
#' @param TypePCA Type of PCA: "PCA" or "SPCA"
#' @param PCA.Files path for PCA file
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param nbclusters number of clusters defined in k-Means
#'
#' @export
map_spectral_species <- function(Input.Image.File, Output.Dir, PCA.Files, TypePCA = "SPCA", nbclusters = 50, nbCPU = 1, MaxRAM = FALSE) {

  # for each image
  Output.Dir2 <- define_output_directory(Output.Dir, Input.Image.File, TypePCA)
  Output.Dir.SS <- define_output_subdir(Output.Dir, Input.Image.File, TypePCA, "SpectralSpecies")
  Output.Dir.PCA <- define_output_subdir(Output.Dir, Input.Image.File, TypePCA, "PCA")
  Spectral.Species.Path <- paste(Output.Dir.SS, "SpectralSpecies", sep = "")
  # if no prior diversity map has been produced --> need PCA file
  if (!file.exists(PCA.Files)) {
    message("")
    message("*********************************************************")
    message("WARNING: PCA file expected but not found")
    print(PCA.Files)
    message("process aborted")
    message("*********************************************************")
    message("")
    stop()
  } else {
    WS_Save <- paste(Output.Dir2, "PCA_Info.RData", sep = "")
    load(file = WS_Save)
    ##            1- Select components used to perform clustering           ##
    PC.Select.Path <- paste(Output.Dir.PCA, "Selected_Components.txt", sep = "")
    if (file.exists(PC.Select.Path)) {
      PC.Select <- read.table(PC.Select.Path)[[1]]
      dataPCA <- PCA.model$dataPCA[, PC.Select]
      if (length(PC.Select) == 1) {
        dataPCA <- matrix(dataPCA, ncol = 1)
      }
      message("Selected components:")
      print(PC.Select)
      message("Please add carriage return after last selected component if not part of the list")
      message("If these do not match with your selection, please correct file following file:")
      print(PC.Select.Path)
    } else {
      print(paste("File named ", PC.Select.Path, "needs to be created first"))
      print("Image processing aborted")
      exit()
    }
    ##    2- PERFORM KMEANS FOR EACH ITERATION & DEFINE SPECTRAL SPECIES    ##
    print("perform k-means clustering for each subset and define centroids")
    # scaling factor subPCA between 0 and 1
    Kmeans.info <- init_kmeans(dataPCA, Pix.Per.Iter, nb_partitions, nbclusters, nbCPU)
    ##              3- ASSIGN SPECTRAL SPECIES TO EACH PIXEL                ##
    print("apply Kmeans to the whole image and determine spectral species")
    apply_kmeans(PCA.Files, PC.Select, ImPathShade, Kmeans.info, Spectral.Species.Path, nbCPU, MaxRAM)
  }
  return()
}

# computes k-means from nb_partitions subsets taken from dataPCA
#
# @param dataPCA initial dataset sampled from PCA image
# @param Pix.Per.Iter number of pixels per iteration
# @param nb_partitions number of k-means then averaged
# @param nbCPU
# @param nbclusters number of clusters used in kmeans
#
# @return list of centroids and parameters needed to center/reduce data
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats kmeans
init_kmeans <- function(dataPCA, Pix.Per.Iter, nb_partitions, nbclusters, nbCPU = 1) {
  m0 <- apply(dataPCA, 2, function(x) min(x))
  M0 <- apply(dataPCA, 2, function(x) max(x))
  d0 <- M0 - m0
  dataPCA <- center_reduce(dataPCA, m0, d0)
  # get the dimensions of the images, and the number of subimages to process
  dataPCA <- split(as.data.frame(dataPCA), rep(1:nb_partitions, each = Pix.Per.Iter))
  # multiprocess kmeans clustering
  plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
  Schedule.Per.Thread <- ceiling(length(dataPCA) / nbCPU)
  res <- future_lapply(dataPCA, FUN = kmeans, centers = nbclusters, iter.max = 50, nstart = 10,
                       algorithm = c("Hartigan-Wong"), future.scheduling = Schedule.Per.Thread)
  plan(sequential)
  Centroids <- list(nb_partitions)
  for (i in (1:nb_partitions)) {
    Centroids[[i]] <- res[[i]]$centers
  }
  my_list <- list("Centroids" = Centroids, "MinVal" = m0, "MaxVal" = M0, "Range" = d0)
  return(my_list)
}

# Applies Kmeans clustering to PCA image
#
# @param PCA.Path path for the PCA image
# @param PC.Select PCs selected from PCA
# @param ImPathShade shade mask
# @param Kmeans.info information about kmeans computed in previous step
# @param nbCPU
# @param MaxRAM
# @param Spectral.Species.Path path for spectral species file to be written
#
# @return
apply_kmeans <- function(PCA.Path, PC.Select, ImPathShade, Kmeans.info, Spectral.Species.Path, nbCPU = 1, MaxRAM = FALSE) {
  ImPathHDR <- get_HDR_name(PCA.Path)
  HDR.PCA <- read_ENVI_header(ImPathHDR)
  PCA.Format <- ENVI_type2bytes(HDR.PCA)
  HDR.Shade <- get_HDR_name(ImPathShade)
  HDR.Shade <- read_ENVI_header(HDR.Shade)
  # prepare for sequential processing: SeqRead.Image informs about byte location to read
  if (MaxRAM == FALSE) {
    MaxRAM <- 0.25
  }
  nbPieces <- split_image(HDR.PCA, MaxRAM)
  SeqRead.PCA <- where_to_read(HDR.PCA, nbPieces)
  SeqRead.Shade <- where_to_read(HDR.Shade, nbPieces)
  # create output file for spectral species assignment
  HDR.SS <- HDR.PCA
  nb_partitions <- length(Kmeans.info$Centroids)
  HDR.SS$bands <- nb_partitions
  HDR.SS$`data type` <- 1
  Iter <- paste('Iter', 1:nb_partitions, collapse = ", ")
  HDR.SS$`band names` <- Iter
  HDR.SS$wavelength <- NULL
  HDR.SS$fwhm <- NULL
  HDR.SS$resolution <- NULL
  HDR.SS$bandwidth <- NULL
  HDR.SS$purpose <- NULL
  HDR.SS$`byte order` <- get_byte_order()
  headerFpath <- paste(Spectral.Species.Path, ".hdr", sep = "")
  write_ENVI_header(HDR.SS, headerFpath)
  # create Spectral species file
  fidSS <- file(
    description = Spectral.Species.Path, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidSS)

  for (i in 1:nbPieces) {
    print(paste("Spectral Species Piece #", i, "/", nbPieces))
    Location.RW <- list()
    Location.RW$nbLines <- SeqRead.PCA$Lines.Per.Chunk[i]
    Location.RW$Byte.Start.PCA <- SeqRead.PCA$ReadByte.Start[i]
    Location.RW$lenBin.PCA <- SeqRead.PCA$ReadByte.End[i] - SeqRead.PCA$ReadByte.Start[i] + 1
    Location.RW$Byte.Start.Shade <- SeqRead.Shade$ReadByte.Start[i]
    Location.RW$lenBin.Shade <- SeqRead.Shade$ReadByte.End[i] - SeqRead.Shade$ReadByte.Start[i] + 1
    Location.RW$Byte.Start.SS <- 1 + (SeqRead.Shade$ReadByte.Start[i] - 1) * nb_partitions
    Location.RW$lenBin.SS <- nb_partitions * (SeqRead.Shade$ReadByte.End[i] - SeqRead.Shade$ReadByte.Start[i]) + 1
    compute_spectral_species(PCA.Path, ImPathShade, Spectral.Species.Path, Location.RW, PC.Select, Kmeans.info, nbCPU)
  }
  return()
}

# this function reads PCA file and defines the spectral species for each pixel
# based on the set of cluster centroids defined for each iteration
# applies kmeans --> closest cluster corresponds to the "spectral species"
#
# @param PCA.Path path for the PCA image
# @param ImPathShade shade mask
# @param Spectral.Species.Path path for spectral species file to be written
# @param Location.RW where to read/write information in binary file
# @param PC.Select PCs selected from PCA
# @param nbCPU
# @param Kmeans.info information about kmeans computed in previous step
#
# @return
#' @importFrom snow splitRows
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
compute_spectral_species <- function(PCA.Path, ImPathShade, Spectral.Species.Path, Location.RW, PC.Select, Kmeans.info, nbCPU = 1) {

  # characteristics of the centroids computed during preprocessing
  nb_partitions <- length(Kmeans.info$Centroids)
  nbCentroids <- nrow(Kmeans.info$Centroids[[1]])
  CentroidsArray <- do.call("rbind", Kmeans.info$Centroids)

  # read shade file and PCA file
  ShadeHDR <- get_HDR_name(ImPathShade)
  HDR.Shade <- read_ENVI_header(ShadeHDR)
  Shade.Format <- ENVI_type2bytes(HDR.Shade)
  ImgFormat <- "Shade"
  Shade.Chunk <- read_image_subset(ImPathShade, HDR.Shade, Location.RW$Byte.Start.Shade, Location.RW$lenBin.Shade, Location.RW$nbLines, Shade.Format, ImgFormat)

  PCA.HDR <- get_HDR_name(PCA.Path)
  HDR.PCA <- read_ENVI_header(PCA.HDR)
  PCA.Format <- ENVI_type2bytes(HDR.PCA)
  # read "unfolded" (2D) PCA image
  ImgFormat <- "2D"
  PCA.Chunk <- read_image_subset(PCA.Path, HDR.PCA, Location.RW$Byte.Start.PCA, Location.RW$lenBin.PCA, Location.RW$nbLines, PCA.Format, ImgFormat)
  PCA.Chunk <- PCA.Chunk[, PC.Select]
  if (length(PC.Select) == 1) {
    PCA.Chunk <- matrix(PCA.Chunk, ncol = 1)
  }
  PCA.Chunk <- center_reduce(PCA.Chunk, Kmeans.info$MinVal, Kmeans.info$Range)
  # eliminate shaded pixels
  keepShade <- which(Shade.Chunk == 1)
  PCA.Chunk <- matrix(PCA.Chunk[keepShade, ], ncol = length(PC.Select))

  # Prepare writing of spectral species distribution file
  SS.HDR <- get_HDR_name(Spectral.Species.Path)
  HDR.SS <- read_ENVI_header(SS.HDR)
  SS.Format <- ENVI_type2bytes(HDR.SS)

  # for each pixel in the subset, compute the nearest cluster for each iteration
  Nearest.Cluster <- matrix(0, nrow = Location.RW$nbLines * HDR.PCA$samples, ncol = nb_partitions)
  # rdist consumes RAM  --> divide data into pieces if too big and multiprocess
  nbSamples.Per.Rdist <- 25000
  if (length(keepShade) > 0) {
    nbSubsets <- ceiling(length(keepShade) / nbSamples.Per.Rdist)
    PCA.Chunk <- splitRows(PCA.Chunk, nbSubsets)

    plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
    Schedule.Per.Thread <- ceiling(nbSubsets / nbCPU)
    res <- future_lapply(PCA.Chunk, FUN = RdistList, CentroidsArray = CentroidsArray, nbClusters = nrow(Kmeans.info$Centroids[[1]]), nb_partitions = nb_partitions, future.scheduling = Schedule.Per.Thread)
    plan(sequential)
    res <- do.call("rbind", res)
    Nearest.Cluster[keepShade, ] <- res
    rm(res)
    gc()
  }
  Nearest.Cluster <- array(Nearest.Cluster, c(Location.RW$nbLines, HDR.PCA$samples, nb_partitions))
  # Write spectral species distribution in the output file
  fidSS <- file(
    description = Spectral.Species.Path, open = "r+b", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  Nearest.Cluster <- aperm(Nearest.Cluster, c(2, 3, 1))
  if (!Location.RW$Byte.Start.SS == 1) {
    seek(fidSS, where = SS.Format$Bytes * (Location.RW$Byte.Start.SS - 1), origin = "start", rw = "write")
  }
  writeBin(c(as.integer(Nearest.Cluster)), fidSS, size = SS.Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
  close(fidSS)
  rm(list = ls())
  gc()
  return()
}

# Compute distance between each pixel of input data and each of the nbClusters x nb_partitions centroids
#
# @param InputData
# @param CentroidsArray
# @param nbClusters
# @param nb_partitions
#
# @return ResDist
#' @importFrom fields rdist
RdistList <- function(InputData, CentroidsArray, nbClusters, nb_partitions) {
  # number of pixels in input data
  nbPixels <- nrow(InputData)
  # compute distance between each pixel and each centroid
  cluster.dist <- rdist(InputData, CentroidsArray)
  # reshape distance into a matrix: all pixels from iteration 1, then all pixels from it2...
  cluster.dist <- matrix(aperm(array(cluster.dist, c(nbPixels, nbClusters, nb_partitions)), c(1, 3, 2)), nrow = nbPixels * nb_partitions)
  ResDist <- matrix(max.col(-cluster.dist), nrow = nbPixels)
  return(ResDist)
}
