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
#' @param Input_Image_File character. Path of the image to be processed
#' @param Output_Dir character. Path for output directory
#' @param TypePCA character. Type of PCA: choose either "PCA" or "SPCA"
#' @param PCA_Files character. Path of the PCA image
#' @param ImPathShade character. Path of the mask corresponding to the image
#' @param Pix_Per_Partition number of pixels for each partition
#' @param nb_partitions number of partition
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param nbclusters number of clusters defined in k-Means
#'
#' @importFrom utils read.table
#' @export
map_spectral_species <- function(Input_Image_File, Output_Dir, PCA_Files,ImPathShade,Pix_Per_Partition,nb_partitions,TypePCA = "SPCA", nbclusters = 50, nbCPU = 1, MaxRAM = FALSE) {

  # for each image
  Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "SpectralSpecies")
  Output_Dir_PCA <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "PCA")
  Spectral_Species_Path <- paste(Output_Dir_SS, "SpectralSpecies", sep = "")
  # if no prior diversity map has been produced --> need PCA file
  if (!file.exists(PCA_Files)) {
    message("")
    message("*********************************************************")
    message("WARNING: PCA file expected but not found")
    print(PCA_Files)
    message("process aborted")
    message("*********************************************************")
    message("")
    stop()
  } else {
    WS_Save <- paste(Output_Dir_PCA, "PCA_Info.RData", sep = "")
    load(file = WS_Save)
    ##            1- Select components used to perform clustering           ##
    PC_Select_Path <- paste(Output_Dir_PCA, "Selected_Components.txt", sep = "")
    if (file.exists(PC_Select_Path)) {
      PC_Select <- read.table(PC_Select_Path)[[1]]
      dataPCA <- PCA_model$dataPCA[, PC_Select]
      if (length(PC_Select) == 1) {
        dataPCA <- matrix(dataPCA, ncol = 1)
      }
      message("Selected components:")
      print(PC_Select)
      message("Please add carriage return after last selected component if not part of the list")
      message("If these do not match with your selection, please correct file following file:")
      print(PC_Select_Path)
    } else {
      print(paste("File named ", PC_Select_Path, "needs to be created first"))
      print("Image processing aborted")
      exit()
    }
    ##    2- PERFORM KMEANS FOR EACH ITERATION & DEFINE SPECTRAL SPECIES    ##
    print("perform k-means clustering for each subset and define centroids")
    # scaling factor subPCA between 0 and 1
    Kmeans_info <- init_kmeans(dataPCA, Pix_Per_Partition, nb_partitions, nbclusters, nbCPU)
    ##              3- ASSIGN SPECTRAL SPECIES TO EACH PIXEL                ##
    print("apply Kmeans to the whole image and determine spectral species")
    apply_kmeans(PCA_Files, PC_Select, ImPathShade, Kmeans_info, Spectral_Species_Path, nbCPU, MaxRAM)
  }
  return()
}

# computes k-means from nb_partitions subsets taken from dataPCA
#
# @param dataPCA initial dataset sampled from PCA image
# @param Pix_Per_Partition number of pixels per iteration
# @param nb_partitions number of k-means then averaged
# @param nbCPU
# @param nbclusters number of clusters used in kmeans
#
# @return list of centroids and parameters needed to center/reduce data
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats kmeans
init_kmeans <- function(dataPCA, Pix_Per_Partition, nb_partitions, nbclusters, nbCPU = 1) {
  m0 <- apply(dataPCA, 2, function(x) min(x))
  M0 <- apply(dataPCA, 2, function(x) max(x))
  d0 <- M0 - m0
  dataPCA <- center_reduce(dataPCA, m0, d0)
  # get the dimensions of the images, and the number of subimages to process
  dataPCA <- split(as.data.frame(dataPCA), rep(1:nb_partitions, each = Pix_Per_Partition))
  # multiprocess kmeans clustering
  plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
  Schedule_Per_Thread <- ceiling(length(dataPCA) / nbCPU)
  res <- future_lapply(dataPCA, FUN = kmeans, centers = nbclusters, iter.max = 50, nstart = 10,
                       algorithm = c("Hartigan-Wong"), future.scheduling = Schedule_Per_Thread)
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
# @param PCA_Path path for the PCA image
# @param PC_Select PCs selected from PCA
# @param ImPathShade shade mask
# @param Kmeans_info information about kmeans computed in previous step
# @param nbCPU
# @param MaxRAM
# @param Spectral_Species_Path path for spectral species file to be written
#
# @return
apply_kmeans <- function(PCA_Path, PC_Select, ImPathShade, Kmeans_info, Spectral_Species_Path, nbCPU = 1, MaxRAM = FALSE) {
  ImPathHDR <- get_HDR_name(PCA_Path)
  HDR_PCA <- read_ENVI_header(ImPathHDR)
  PCA_Format <- ENVI_type2bytes(HDR_PCA)
  HDR_Shade <- get_HDR_name(ImPathShade)
  HDR_Shade <- read_ENVI_header(HDR_Shade)
  # prepare for sequential processing: SeqRead_Image informs about byte location to read
  if (MaxRAM == FALSE) {
    MaxRAM <- 0.25
  }
  nbPieces <- split_image(HDR_PCA, MaxRAM)
  SeqRead_PCA <- where_to_read(HDR_PCA, nbPieces)
  SeqRead_Shade <- where_to_read(HDR_Shade, nbPieces)
  # create output file for spectral species assignment
  HDR_SS <- HDR_PCA
  nb_partitions <- length(Kmeans_info$Centroids)
  HDR_SS$bands <- nb_partitions
  HDR_SS$`data type` <- 1
  Iter <- paste('Iter', 1:nb_partitions, collapse = ", ")
  HDR_SS$`band names` <- Iter
  HDR_SS$wavelength <- NULL
  HDR_SS$fwhm <- NULL
  HDR_SS$resolution <- NULL
  HDR_SS$bandwidth <- NULL
  HDR_SS$purpose <- NULL
  HDR_SS$`byte order` <- get_byte_order()
  headerFpath <- paste(Spectral_Species_Path, ".hdr", sep = "")
  write_ENVI_header(HDR_SS, headerFpath)
  # create Spectral species file
  fidSS <- file(
    description = Spectral_Species_Path, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidSS)

  for (i in 1:nbPieces) {
    print(paste("Spectral Species Piece #", i, "/", nbPieces))
    Location_RW <- list()
    Location_RW$nbLines <- SeqRead_PCA$Lines_Per_Chunk[i]
    Location_RW$Byte_Start_PCA <- SeqRead_PCA$ReadByte_Start[i]
    Location_RW$lenBin_PCA <- SeqRead_PCA$ReadByte_End[i] - SeqRead_PCA$ReadByte_Start[i] + 1
    Location_RW$Byte_Start_Shade <- SeqRead_Shade$ReadByte_Start[i]
    Location_RW$lenBin_Shade <- SeqRead_Shade$ReadByte_End[i] - SeqRead_Shade$ReadByte_Start[i] + 1
    Location_RW$Byte_Start_SS <- 1 + (SeqRead_Shade$ReadByte_Start[i] - 1) * nb_partitions
    Location_RW$lenBin_SS <- nb_partitions * (SeqRead_Shade$ReadByte_End[i] - SeqRead_Shade$ReadByte_Start[i]) + 1
    compute_spectral_species(PCA_Path, ImPathShade, Spectral_Species_Path, Location_RW, PC_Select, Kmeans_info, nbCPU)
  }
  return()
}

# this function reads PCA file and defines the spectral species for each pixel
# based on the set of cluster centroids defined for each iteration
# applies kmeans --> closest cluster corresponds to the "spectral species"
#
# @param PCA_Path path for the PCA image
# @param ImPathShade shade mask
# @param Spectral_Species_Path path for spectral species file to be written
# @param Location_RW where to read/write information in binary file
# @param PC_Select PCs selected from PCA
# @param nbCPU
# @param Kmeans_info information about kmeans computed in previous step
#
# @return
#' @importFrom snow splitRows
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
compute_spectral_species <- function(PCA_Path, ImPathShade, Spectral_Species_Path, Location_RW, PC_Select, Kmeans_info, nbCPU = 1) {

  # characteristics of the centroids computed during preprocessing
  nb_partitions <- length(Kmeans_info$Centroids)
  nbCentroids <- nrow(Kmeans_info$Centroids[[1]])
  CentroidsArray <- do.call("rbind", Kmeans_info$Centroids)

  # read shade file and PCA file
  ShadeHDR <- get_HDR_name(ImPathShade)
  HDR_Shade <- read_ENVI_header(ShadeHDR)
  Shade.Format <- ENVI_type2bytes(HDR_Shade)
  ImgFormat <- "Shade"
  Shade_Chunk <- read_image_subset(ImPathShade, HDR_Shade, Location_RW$Byte_Start_Shade, Location_RW$lenBin_Shade, Location_RW$nbLines, Shade.Format, ImgFormat)

  PCA_HDR <- get_HDR_name(PCA_Path)
  HDR_PCA <- read_ENVI_header(PCA_HDR)
  PCA_Format <- ENVI_type2bytes(HDR_PCA)
  # read "unfolded" (2D) PCA image
  ImgFormat <- "2D"
  PCA_Chunk <- read_image_subset(PCA_Path, HDR_PCA, Location_RW$Byte_Start_PCA, Location_RW$lenBin_PCA, Location_RW$nbLines, PCA_Format, ImgFormat)
  PCA_Chunk <- PCA_Chunk[, PC_Select]
  if (length(PC_Select) == 1) {
    PCA_Chunk <- matrix(PCA_Chunk, ncol = 1)
  }
  PCA_Chunk <- center_reduce(PCA_Chunk, Kmeans_info$MinVal, Kmeans_info$Range)
  # eliminate shaded pixels
  keepShade <- which(Shade_Chunk == 1)
  PCA_Chunk <- matrix(PCA_Chunk[keepShade, ], ncol = length(PC_Select))

  # Prepare writing of spectral species distribution file
  SS_HDR <- get_HDR_name(Spectral_Species_Path)
  HDR_SS <- read_ENVI_header(SS_HDR)
  SS_Format <- ENVI_type2bytes(HDR_SS)

  # for each pixel in the subset, compute the nearest cluster for each iteration
  Nearest_Cluster <- matrix(0, nrow = Location_RW$nbLines * HDR_PCA$samples, ncol = nb_partitions)
  # rdist consumes RAM  --> divide data into pieces if too big and multiprocess
  nbSamples_Per_Rdist <- 25000
  if (length(keepShade) > 0) {
    nbSubsets <- ceiling(length(keepShade) / nbSamples_Per_Rdist)
    PCA_Chunk <- splitRows(PCA_Chunk, nbSubsets)

    plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
    Schedule_Per_Thread <- ceiling(nbSubsets / nbCPU)
    res <- future_lapply(PCA_Chunk, FUN = RdistList, CentroidsArray = CentroidsArray, nbClusters = nrow(Kmeans_info$Centroids[[1]]), nb_partitions = nb_partitions, future.scheduling = Schedule_Per_Thread)
    plan(sequential)
    res <- do.call("rbind", res)
    Nearest_Cluster[keepShade, ] <- res
    rm(res)
    gc()
  }
  Nearest_Cluster <- array(Nearest_Cluster, c(Location_RW$nbLines, HDR_PCA$samples, nb_partitions))
  # Write spectral species distribution in the output file
  fidSS <- file(
    description = Spectral_Species_Path, open = "r+b", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  Nearest_Cluster <- aperm(Nearest_Cluster, c(2, 3, 1))
  if (!Location_RW$Byte_Start_SS == 1) {
    seek(fidSS, where = SS_Format$Bytes * (Location_RW$Byte_Start_SS - 1), origin = "start", rw = "write")
  }
  writeBin(c(as.integer(Nearest_Cluster)), fidSS, size = SS_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
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
  cluster_dist <- rdist(InputData, CentroidsArray)
  # reshape distance into a matrix: all pixels from iteration 1, then all pixels from it2...
  cluster_dist <- matrix(aperm(array(cluster_dist, c(nbPixels, nbClusters, nb_partitions)), c(1, 3, 2)), nrow = nbPixels * nb_partitions)
  ResDist <- matrix(max.col(-cluster_dist), nrow = nbPixels)
  return(ResDist)
}
