# ==============================================================================
# biodivMapR
# Lib_MapSpectralSpecies.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library applies clustering on a selection of components stored in a PCA
# file previously created ("Perform_PCA.R") and produces corresponding spectral
# species
# ==============================================================================

#' maps spectral species based on PCA file computed previously
#'
#' @param Input_Image_File character. Path of the image to be processed
#' @param Input_Mask_File character. Path of the mask corresponding to the image
#' @param Output_Dir character. Path for output directory
#' @param SpectralSpace_Output list. list of variables produced from function perform_PCA
#' @param nbclusters numeric. number of clusters defined in k-Means
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param Kmeans_Only boolean. set to TRUE if computation of kmeans without production of spectral species map
#' @param SelectedPCs numeric. Define PCs to be selected. Set to FALSE if you want to use the "Selected_Components.txt" file
#' @param SpectralFilter list. information about spectral band location
#' (central wavelength), bands to keep...
#'
#' @return Kmeans_info
#' @importFrom utils read.table
#' @export
map_spectral_species <- function(Input_Image_File, Input_Mask_File, Output_Dir,
                                 SpectralSpace_Output, nbclusters = 50,
                                 nbCPU = 1, MaxRAM = 0.25,
                                 Kmeans_Only = FALSE, SelectedPCs = FALSE,
                                 SpectralFilter = NULL) {

  Kmeans_info <- NULL
  # if no prior diversity map has been produced --> need PCA file
  if (!file.exists(SpectralSpace_Output$PCA_Files)) {
    error_no_PCA_file(SpectralSpace_Output$PCA_Files)
    stop()
  }

  # define directories
  Output_Dir_SS <- define_output_subdir(Output_Dir, Input_Image_File, SpectralSpace_Output$TypePCA, "SpectralSpecies")
  Output_Dir_PCA <- define_output_subdir(Output_Dir, Input_Image_File, SpectralSpace_Output$TypePCA, "PCA")
  Spectral_Species_Path <-  file.path(Output_Dir_SS, "SpectralSpecies")

  # 1- Select components used to perform clustering
  if (SelectedPCs == FALSE){
    PC_Select_Path <- file.path(Output_Dir_PCA, "Selected_Components.txt")
  } else {
    PC_Select_Path = 'NoFile'
  }
  if (file.exists(PC_Select_Path)) {
    PC_Select <- utils::read.table(PC_Select_Path)[[1]]
  } else if (!SelectedPCs == FALSE){
    PC_Select <- SelectedPCs
  } else {
    error_PC_sel(Output_Dir_PCA)
    stop()
  }
  message("Selected components:")
  print(PC_Select)

  # 2- sample data from PCA image
  ImNames <- list(Input_Image = Input_Image_File,
                  Mask_list = Input_Mask_File)
  if (is.null(SpectralSpace_Output$nb_partitions)){
    nb_partitions <- 20
  } else {
    nb_partitions <- SpectralSpace_Output$nb_partitions
  }
  Pix_Per_Partition <- define_pixels_per_iter(ImNames, nb_partitions = nb_partitions)

  ImPathHDR <- get_HDR_name(SpectralSpace_Output$PCA_Files)
  HDR <- read_ENVI_header(ImPathHDR)
  Subset <- get_random_subset_from_image(ImPath = SpectralSpace_Output$PCA_Files,
                                         MaskPath = Input_Mask_File,
                                         nb_partitions = nb_partitions,
                                         Pix_Per_Partition = Pix_Per_Partition,
                                         kernel = NULL,MaxRAM = MaxRAM)
  SubsetInit <- Subset
  dataPCA <- Subset$DataSubset[, PC_Select]
  if (length(PC_Select) == 1) {
    dataPCA <- matrix(dataPCA, ncol = 1)
  }

  # 3- PERFORM KMEANS FOR EACH ITERATION & DEFINE SPECTRAL SPECIES
  print("perform k-means clustering for each subset and define centroids")
  Kmeans_info <- init_kmeans(dataPCA = dataPCA,
                             nb_partitions = nb_partitions,
                             nbclusters = nbclusters,
                             nbCPU = nbCPU)
  Kmeans_info$SpectralSpecies <- Spectral_Species_Path

  if (Kmeans_info$Error==FALSE){
    if (Kmeans_Only==FALSE){
      ##    3- ASSIGN SPECTRAL SPECIES TO EACH PIXEL
      apply_kmeans(PCA_Path = SpectralSpace_Output$PCA_Files,
                   PC_Select = PC_Select,
                   Input_Mask_File = Input_Mask_File,
                   Kmeans_info = Kmeans_info,
                   Spectral_Species_Path = Spectral_Species_Path,
                   nbCPU = nbCPU, MaxRAM = MaxRAM)
    } else {
      print("'Kmeans_Only' was set to TRUE: kmeans was not applied on the full image")
      print("Please set 'Kmeans_Only' to FALSE if you want to produce spectral species map")
    }
    # save kmeans info into binary variable
    Kmeans_Path <- file.path(Output_Dir_PCA, "Kmeans_Info.RData")
    save(Kmeans_info, file = Kmeans_Path)
  } else {
    ##    produce error report
    # create directory where error should be stored
    Output_Dir_Error <- define_output_subdir(Output_Dir, Input_Image_File, SpectralSpace_Output$TypePCA, "ErrorReport")
    # identify which samples cause problems
    LocError <- unique(c(which(!is.finite(Kmeans_info$MinVal)),which(!is.finite(Kmeans_info$MaxVal))))
    ValError <- which(!is.finite(dataPCA[,LocError[1]]))
    # Get the original data corresponding to the first sample
    DataError <- SubsetInit$DataSubset[ValError,]
    DataErrorCR <- Subset$DataSubset[ValError,]
    CoordinatesError <- SubsetInit$coordPix[ValError,]
    # save these in a RData file
    FileError <- file.path(Output_Dir_Error,'ErrorPixels.RData')
    ErrorReport <- list('CoordinatesError' = CoordinatesError,'DataError' = DataError,
                        'DataError_afterCR' = DataErrorCR, 'SpectralFilter'=SpectralFilter)
    save(ErrorReport, file = FileError)
    message("")
    message("*********************************************************")
    message("       An error report directory has been produced.      ")
    message("Please check information about data causing errors here: ")
    print(FileError)
    message("               The process will now stop                 ")
    message("*********************************************************")
    message("")
    stop()
  }
  return(Kmeans_info)
}

#' computes k-means from nb_partitions subsets taken from dataPCA
#'
#' @param dataPCA numeric. initial dataset sampled from PCA image
#' @param nb_partitions numeric. number of k-means then averaged
#' @param nbCPU numeric. Number of CPUs available
#' @param nbclusters numeric.number of clusters used in kmeans
#'
#' @return list of centroids and parameters needed to center/reduce data
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats kmeans
#' @importFrom snow splitRows
#'
#'
#' @export

init_kmeans <- function(dataPCA, nb_partitions, nbclusters, nbCPU = 1) {

  # define boundaries defining outliers based on IQR
  m0 <- M0 <- c()
  for (i in 1:ncol(dataPCA)){
    IQR <- IQR_outliers(dataPCA[,i], weightIRQ = 2)
    m0 <- c(m0,IQR[1])
    M0 <- c(M0,IQR[2])
  }
  d0 <- M0 - m0

  # m0 <- apply(dataPCA, 2, function(x) min(x))
  # M0 <- apply(dataPCA, 2, function(x) max(x))
  # d0 <- M0 - m0
  if (length(which(is.na(m0)))>0 | length(which(is.na(M0)))>0 | length(which(is.infinite(m0)))>0 | length(which(is.infinite(M0)))>0){
    error_input()
    my_list <- list("Centroids" = NULL, "MinVal" = m0, "MaxVal" = M0, "Range" = d0, "Error" = TRUE)
    return(my_list)
  } else {
    dataPCA <- center_reduce(dataPCA, m0, d0)
    dataPCA <- snow::splitRows(x = dataPCA, ncl = nb_partitions)
    if (nbCPU>1){
      plan(multisession, workers = nbCPU) ## Parallelize using four cores
      handlers(global = TRUE)
      handlers("cli")
      with_progress({
        p <- progressr::progressor(steps = nb_partitions)
        res <- future_lapply(X = dataPCA,
                             FUN = kmeans_progressr,
                             centers = nbclusters,
                             iter.max = 50, nstart = 10,
                             algorithm = c("Hartigan-Wong"), p = p)
      })
      plan(sequential)
    } else {
      handlers(global = TRUE)
      handlers("cli")
      with_progress({
        p <- progressr::progressor(steps = nb_partitions)
        res <- lapply(X = dataPCA,
                      FUN = kmeans_progressr,
                      centers = nbclusters, iter.max = 50, nstart = 10,
                      algorithm = c("Hartigan-Wong"), p = p)
      })
    }
    Centroids <- lapply(res,'[[',2)
    my_list <- list("Centroids" = Centroids, "MinVal" = m0, "MaxVal" = M0, "Range" = d0, "Error" = FALSE)
    return(my_list)
  }
}

#' applies results of ordination to full image based on nearest neighbors
#
#' @param x numeric. subset of spectral species distribution file
#' @param centers numeric. Samples selected during ordination
#' @param iter.max numeric. ordination of dissimilarity matrix for a selection of spatial units
#' @param nstart numeric. Number of partitions (repetitions) to be computed then averaged.
#' @param algorithm character. Name of teh algorithm
#' @param p function.
#
#' @return results of kmeans
#' @importFrom stats kmeans
#' @export

kmeans_progressr <- function(x, centers, iter.max, nstart,
                             algorithm = c("Hartigan-Wong"), p = NULL){
  res <- kmeans(x = x, centers = centers, iter.max = iter.max, nstart = nstart)
  if (!is.null(p)){p()}
  return(res)
}

#' Applies Kmeans clustering to PCA image and writes spectral species map
#'
#' @param PCA_Path path for the PCA image
#' @param PC_Select PCs selected from PCA
#' @param Input_Mask_File Path for the mask
#' @param Kmeans_info information about kmeans computed in previous step
#' @param Spectral_Species_Path path for spectral species file to be written
#' @param nbCPU numeric. number of CPU to work with in multisession task
#' @param MaxRAM numeric. size of image pieces to be read at once
#'
#' @return None
#' @importFrom progress progress_bar
#' @export

apply_kmeans <- function(PCA_Path, PC_Select, Input_Mask_File, Kmeans_info,
                         Spectral_Species_Path, nbCPU = 1, MaxRAM = 0.25) {

  nb_partitions <- length(Kmeans_info$Centroids)

  HDR_PCA <- read_ENVI_header(get_HDR_name(PCA_Path))
  PCA_Format <- ENVI_type2bytes(HDR_PCA)
  HDR_Shade <- read_ENVI_header(get_HDR_name(Input_Mask_File))
  # prepare for sequential processing: SeqRead_Image informs about byte location to read
  nbPieces <- split_image(HDR_PCA, MaxRAM)
  SeqRead_PCA <- where_to_read(HDR_PCA, nbPieces)
  SeqRead_Shade <- where_to_read(HDR_Shade, nbPieces)

  Location_RW <- list()
  for (i in 1:nbPieces) {
    Location_RW[[i]] <- list()
    Location_RW[[i]]$nbLines <- SeqRead_PCA$Lines_Per_Chunk[i]
    # PCA file
    Location_RW[[i]]$Byte_Start_PCA <- SeqRead_PCA$ReadByte_Start[i]
    Location_RW[[i]]$lenBin_PCA <- SeqRead_PCA$ReadByte_End[i] - SeqRead_PCA$ReadByte_Start[i] + 1
    # shade file
    Location_RW[[i]]$Byte_Start_Shade <- SeqRead_Shade$ReadByte_Start[i]
    Location_RW[[i]]$lenBin_Shade <- SeqRead_Shade$ReadByte_End[i] - SeqRead_Shade$ReadByte_Start[i] + 1
    # spectral species file
    Location_RW[[i]]$Byte_Start_SS <- 1 + (SeqRead_Shade$ReadByte_Start[i] - 1) * nb_partitions
    Location_RW[[i]]$lenBin_SS <- nb_partitions * (SeqRead_Shade$ReadByte_End[i] - SeqRead_Shade$ReadByte_Start[i]) + 1
  }

  # create output file for spectral species assignment
  create_HDR_SS(HDR_tempate = HDR_PCA,
                nbBands = nb_partitions,
                Spectral_Species_Path)
  # create Spectral species file
  fidSS <- file(
    description = Spectral_Species_Path, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidSS)

  # for each piece of image
  print(paste('Apply Kmeans to the full raster:',nbPieces,'chunks distributed on',nbCPU,'CPU'))
  # pb <- progress_bar$new(
  #   format = paste('Write spectral species file [:bar] :percent in :elapsedfull',sep = ''),
  #   total = nbPieces, clear = FALSE, width= 100)
  # lapply(X = Location_RW, FUN = compute_spectral_species,
  #        PCA_Path = PCA_Path,
  #        Input_Mask_File = Input_Mask_File,
  #        Spectral_Species_Path = Spectral_Species_Path,
  #        PC_Select = PC_Select,
  #        Kmeans_info = Kmeans_info,
  #        nbCPU = nbCPU)

  for (i in 1:nbPieces) {
    message(paste('Computing spectral species for image subset #',i,' / ',nbPieces))
    compute_spectral_species(PCA_Path = PCA_Path, Input_Mask_File = Input_Mask_File,
                             Spectral_Species_Path = Spectral_Species_Path, Location_RW = Location_RW[[i]],
                             PC_Select = PC_Select, Kmeans_info = Kmeans_info, nbCPU = nbCPU)
  }
  return(invisible())
}

#' this function reads PCA file and defines the spectral species for each pixel
#' based on the set of cluster centroids defined for each iteration
#' applies kmeans --> closest cluster corresponds to the "spectral species"
#'
#' @param PCA_Path character. path for the PCA image
#' @param Input_Mask_File character. Path for the mask
#' @param Spectral_Species_Path character. path for spectral species file to be written
#' @param Location_RW numeric. where to read/write information in binary file
#' @param PC_Select numeric. PCs selected from PCA
#' @param Kmeans_info list. information about kmeans computed in previous step
#' @param nbCPU numeric. number of CPUs available
#'
#' @return None
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom snow splitRows
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @export

compute_spectral_species <- function(PCA_Path, Input_Mask_File, Spectral_Species_Path,
                                     Location_RW, PC_Select, Kmeans_info, nbCPU = 1) {

  # characteristics of the centroids computed during preprocessing
  nb_partitions <- length(Kmeans_info$Centroids)
  nbCentroids <- nrow(Kmeans_info$Centroids[[1]])
  CentroidsArray <- do.call("rbind", Kmeans_info$Centroids)

  # read shade file and PCA file
  ShadeHDR <- get_HDR_name(Input_Mask_File)
  HDR_Shade <- read_ENVI_header(ShadeHDR)
  Shade.Format <- ENVI_type2bytes(HDR_Shade)
  ImgFormat <- "Shade"
  Shade_Chunk <- read_BIL_image_subset(Input_Mask_File, HDR_Shade, Location_RW$Byte_Start_Shade,
                                       Location_RW$lenBin_Shade, Location_RW$nbLines, Shade.Format, ImgFormat)

  PCA_HDR <- get_HDR_name(PCA_Path)
  HDR_PCA <- read_ENVI_header(PCA_HDR)
  PCA_Format <- ENVI_type2bytes(HDR_PCA)
  # read "unfolded" (2D) PCA image
  ImgFormat <- "2D"
  PCA_Chunk <- read_BIL_image_subset(PCA_Path, HDR_PCA, Location_RW$Byte_Start_PCA, Location_RW$lenBin_PCA,
                                     Location_RW$nbLines, PCA_Format, ImgFormat)
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
  # rdist consumes RAM  --> divide data into pieces if too big and multisession
  nbSamples_Per_Rdist <- 25000
  if (length(keepShade) > 0) {
    nbSubsets <- ceiling(length(keepShade) / nbSamples_Per_Rdist)
    PCA_Chunk <- snow::splitRows(PCA_Chunk, nbSubsets)
    if (nbCPU>1){
      plan(multisession, workers = nbCPU)
      handlers(global = TRUE)
      handlers("cli")
      with_progress({
        p <- progressr::progressor(steps = nbSubsets)
        res <- future_lapply(PCA_Chunk,
                             FUN = RdistList,
                             CentroidsArray = CentroidsArray,
                             nbClusters = nrow(Kmeans_info$Centroids[[1]]),
                             nb_partitions = nb_partitions, p = p)
      })
      plan(sequential)
    } else {
      res <- lapply(PCA_Chunk,
                    FUN = RdistList,
                    CentroidsArray = CentroidsArray,
                    nbClusters = nrow(Kmeans_info$Centroids[[1]]),
                    nb_partitions = nb_partitions, p = NULL)
    }
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
  return(invisible())
}

#' compute spectral species for a subset of pixels provided in a list, each element
#' corresponding to a polygon
#'
#' @param subset_Raster numeric. Subset of a raster file on which computation of spectral species sould be performed
#' @param List_FieldPlot list. list of information from same file as subset, corresponding to field plots
#' @param nb_partitions numeric. number of repetitions of kmeans
#' @param Pix_Per_Partition numeric. number of pixels per partition
#' @param nbclusters numeric. number of clusters / spectral species
#' @param PC_Select numeric. selection of components from subset_Raster and List_FieldPlot on which spctral species are computed
#' @param nbCPU numeric. number of CPU on which kmeans is computed
#'
#' @return list. vector_coordinates and vector_ID for each element in the vector file
#' @export

compute_spectral_species_FieldPlots <- function(subset_Raster, List_FieldPlot, nb_partitions,
                                                Pix_Per_Partition, nbclusters, PC_Select=NULL, nbCPU=1){
  # COMPUTE KMEANS
  nbFieldPlots <- length(List_FieldPlot)
  if (!is.null(PC_Select)){
    subset_Raster <- matrix(subset_Raster[,PC_Select],ncol = length(PC_Select))
    for (ll in 1:nbFieldPlots){
      List_FieldPlot[[ll]] <- matrix(List_FieldPlot[[ll]][, PC_Select],ncol = length(PC_Select))
    }
  }
  Kmeans_info <- init_kmeans(subset_Raster, nb_partitions, nbclusters, nbCPU)
  # APPLY KMEANS ON THE FIELD DATA
  Nearest_Cluster <- list()
  # prepare to save spectral species for each plot, for this band combination
  SpectralSpecies_Plots <- list()
  for (ip in 1:nbFieldPlots){
    # if only one polygon in the shapefile and if the polyon is not included in the Raster_SpectralSpecies
    if (length(List_FieldPlot[[ip]])>0){
      PCA_plot <- center_reduce(List_FieldPlot[[ip]], Kmeans_info$MinVal, Kmeans_info$Range)
      # for each pixel in the subset, compute the nearest cluster for each iteration
      Nearest_Cluster[[ip]] <- matrix(0, nrow = nrow(PCA_plot), ncol = nb_partitions)
      CentroidsArray <- do.call("rbind", Kmeans_info$Centroids)
      Nearest_Cluster[[ip]] <- RdistList(PCA_plot, CentroidsArray, nbclusters, nb_partitions)
      SpectralSpecies_Plots[[ip]] <- Nearest_Cluster[[ip]]
    }
  }
  return(SpectralSpecies_Plots)
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
RdistList <- function(InputData, CentroidsArray, nbClusters, nb_partitions, p = NULL) {
  # number of pixels in input data
  nbPixels <- nrow(InputData)
  # compute distance between each pixel and each centroid
  cluster_dist <- rdist(InputData, CentroidsArray)
  # reshape distance into a matrix: all pixels from iteration 1, then all pixels from it2...
  cluster_dist <- matrix(aperm(array(cluster_dist, c(nbPixels, nbClusters, nb_partitions)), c(1, 3, 2)), nrow = nbPixels * nb_partitions)
  ResDist <- matrix(max.col(-cluster_dist), nrow = nbPixels)
  rm(cluster_dist)
  gc()
  if (!is.null(p)){p()}
  return(ResDist)
}

# create HDR for spectral species file based on the HDR from PCA file
create_HDR_SS <- function(HDR_tempate, nbBands, Spectral_Species_Path){
  HDR_SS <- HDR_tempate
  HDR_SS$bands <- nbBands
  HDR_SS$`data type` <- 1
  HDR_SS$`band names` <- paste('Iter', 1:nbBands, collapse = ", ")
  HDR_SS$wavelength <- HDR_SS$fwhm <- HDR_SS$resolution <- HDR_SS$bandwidth <- NULL
  HDR_SS$`default bands` <- HDR_SS$`wavelength units` <- HDR_SS$`z plot titles` <- NULL
  HDR_SS$purpose <- HDR_SS$`data gain values` <- HDR_SS$`default stretch` <- NULL
  HDR_SS$interleave <- 'BIL'
  HDR_SS$`byte order` <- get_byte_order()
  headerFpath <- paste(Spectral_Species_Path, ".hdr", sep = "")
  write_ENVI_header(HDR_SS, headerFpath)
  return()
}

# prints an error message if problems occur
error_input <- function() {
  message("")
  message("*********************************************************")
  message("WARNING: the processing resulted in NA or infinite values")
  message("     This may be due to noisy spectral domains or        ")
  message(" individual pixels showing Inf or Na values in input data")
  message("               Please check input data                   ")
  message("                                                         ")
  message("   if nothing wrong identified with input data, please   ")
  message("   submit a bug report, reproduce the bug with reduced   ")
  message("     dataset and contact the authors of the package      ")
  message("                   process aborted                       ")
  message("*********************************************************")
  message("")
  return()
}


# prints an error message if no PCA file is found
error_no_PCA_file <- function(PCA_Files){
  message("")
  message("*********************************************************")
  message("WARNING: This file required to compute spectral species is missing")
  print(PCA_Files)
  message("process aborted")
  message("*********************************************************")
  message("")
  return()
}

# prints an error message if no PCA file is found
error_PC_sel <- function(Output_Dir_PCA){
  print("PC SELECTION MUST BE PERFORMED FIRST")
  print("Please identify selected components either in this file:")
  PC_Select_Path <- file.path(Output_Dir_PCA, "Selected_Components.txt")
  print(PC_Select_Path)
  print("or in the 'SelectedPCs' variable of map_spectral_species")
  print("Image processing aborted")
}




