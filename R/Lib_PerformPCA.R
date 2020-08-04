# ==============================================================================
# biodivMapR
# Lib_PerformPCA.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library is used to perform PCA on raster prior to diversity mapping
# ==============================================================================

#' Performs PCA for all images and create PCA file with either all or a selection of PCs
#'
#' @param Input_Image_File character. Path of the image to be processed
#' @param Input_Mask_File character. Path of the mask corresponding to the image
#' @param Output_Dir character. Path for output directory
#' @param Continuum_Removal boolean. Set to TRUE if continuum removal should be applied
#' @param TypePCA character. Type of PCA: choose either "PCA" or "SPCA"
#' @param NbPCs_To_Keep numeric. number of components to ke saved in the PCA file. default = 30 if set to FALSE (or nb PC if <30)
#' @param FilterPCA boolean. Set to TRUE if 2nd filtering based on PCA is required
#' @param Excluded_WL  numeric. Water Vapor Absorption domains (in nanometers, min and max WL). Can also be used to exclude spectific domains. dims = N x 2 (N = number of domains to be eliminated)
#' @param nb_partitions numeric. Number of repetitions to estimate diversity from the raster (averaging repetitions).
#' @param nbCPU numeric. Number fo CPUs in use.
#' @param MaxRAM numeric. Maximum size of chunk in GB to limit RAM allocation when reading image file.
#'
#' @return list of paths corresponding to resulting PCA files
#' @export
perform_PCA  <- function(Input_Image_File, Input_Mask_File, Output_Dir, Continuum_Removal=TRUE, TypePCA="SPCA",
                         NbPCs_To_Keep=30,FilterPCA = FALSE, Excluded_WL = FALSE, nb_partitions = 20,
                         nbCPU = 1, MaxRAM = 0.25) {
  # check if format of raster data is as expected
  check_data(Input_Image_File)
  if (!Input_Mask_File==FALSE){
    check_data(Input_Mask_File,Mask = TRUE)
  }
  # define the path corresponding to image, mask and output directory
  ImNames <- list()
  ImNames$Input.Image <- Input_Image_File
  ImNames$Mask_list <- Input_Mask_File
  Output_Dir_Full <- define_output_directory(Output_Dir, Input_Image_File, TypePCA)
  # Identify water vapor absorption bands in image and possibly other spectral domains to discard
  SpectralFilter <- exclude_spectral_domains(Input_Image_File, Excluded_WL = Excluded_WL)
  # Extract valid data subset and check validity
  print("Extract pixels from the images to perform PCA on a subset")
  # define number of pixels to be extracted from the image for each iteration
  Pix_Per_Partition <- define_pixels_per_iter(ImNames, nb_partitions = nb_partitions)
  nb_Pix_To_Sample <- nb_partitions * Pix_Per_Partition
  ImPathHDR <- get_HDR_name(Input_Image_File)
  HDR <- read_ENVI_header(ImPathHDR)
  # extract a random selection of pixels from image
  if (TypePCA=='MNF'){
    FilterPCA <- FALSE
    kernel = matrix(0, 3, 3)
    kernel[c(5, 6, 8)]=c(1, -1/2, -1/2)
    Subset <- get_random_subset_from_image(ImPath = Input_Image_File,
                                           MaskPath = Input_Mask_File, nb_partitions = nb_partitions,
                                           Pix_Per_Partition = Pix_Per_Partition, kernel = kernel, MaxRAM = MaxRAM)
  } else {
    Subset <- get_random_subset_from_image(ImPath = Input_Image_File,
                                           MaskPath = Input_Mask_File, nb_partitions = nb_partitions,
                                           Pix_Per_Partition = Pix_Per_Partition, kernel = NULL, MaxRAM = MaxRAM)
  }
  # if needed, apply continuum removal
  if (Continuum_Removal == TRUE) {
    Subset$DataSubset <- apply_continuum_removal(Subset$DataSubset, SpectralFilter, nbCPU = nbCPU)
  } else {
    if (!length(SpectralFilter$WaterVapor) == 0) {
      Subset$DataSubset <- Subset$DataSubset[, -SpectralFilter$WaterVapor]
    }
  }
  # if number of pixels available inferior number initial sample size
  if (Subset$nbPix2Sample < nb_Pix_To_Sample) {
    nb_Pix_To_Sample <- Subset$nbPix2Sample
    nb_partitions <- ceiling(nb_Pix_To_Sample / Pix_Per_Partition)
    Pix_Per_Partition <- floor(nb_Pix_To_Sample / nb_partitions)
    nb_Pix_To_Sample <- nb_partitions * Pix_Per_Partition
  }
  DataSubset <- Subset$DataSubset
  # clean reflectance data from inf and constant values
  CleanData <- rm_invariant_bands(DataSubset, SpectralFilter)
  DataSubset <- CleanData$DataMatrix
  SpectralFilter <- CleanData$Spectral

  # Compute PCA #1 on DataSubset
  print("perform PCA#1 on the subset image")
  if (TypePCA == "PCA" | TypePCA == "SPCA") {
    PCA_model <- pca(DataSubset, TypePCA)
  # } else if (TypePCA == "NLPCA") {
  #   print("performing NL-PCA with autoencoder")
  #   print("Make sure you properly installed and defined python environment if using this functionality")
  #   tic()
  #   PCA_model <- nlpca(DataSubset)
  #   toc()
  } else if(TypePCA=="MNF"){
    PCA_model <- mnf(DataSubset, Subset$coordPix)
  }

  # if PCA based filtering:
  if (FilterPCA == TRUE) {
    # Perform PCA-based pixels filtering
    # the shade mask helps discarding most unwanted pixels: Shade, clouds, soil,
    # water...). However some unwanted pixels remain. Here we hypothese that
    # such pixels which do not correspond to vegetation will take extreme values
    # after PCA transformation.
    # In order to exclude these pixels, we compute mean and SD for the 3 first
    # components and exclude all pixels showing values ouside "mean+-3SD" range
    print("perform 2nd filtering: Exclude extreme PCA values")
    if (dim(PCA_model$x)[2] > 5) {
      PCsel <- 1:5
    } else {
      PCsel <- 1:dim(PCA_model$x)[2]
    }
    Shade_Update <- paste(Output_Dir_Full, "ShadeMask_Update_PCA", sep = "")
    Input_Mask_File <- filter_PCA(Input_Image_File, HDR, Input_Mask_File, Shade_Update,
                                  Spectral = SpectralFilter,Continuum_Removal, PCA_model,
                                  PCsel, TypePCA,nbCPU = nbCPU, MaxRAM = MaxRAM)
    ## Compute PCA 2 based on the updated shade mask ##
    # extract a random selection of pixels from image
    Subset <- get_random_subset_from_image(ImPath = Input_Image_File, MaskPath = Input_Mask_File,
                                           nb_partitions = nb_partitions, Pix_Per_Partition = Pix_Per_Partition,
                                           kernel = NULL, MaxRAM = MaxRAM)
    # if needed, apply continuum removal
    if (Continuum_Removal == TRUE) {
      Subset$DataSubset <- apply_continuum_removal(Subset$DataSubset, SpectralFilter, nbCPU = nbCPU)
    } else {
      if (!length(SpectralFilter$WaterVapor) == 0) {
        Subset$DataSubset <- Subset$DataSubset[, -SpectralFilter$WaterVapor]
      }
    }
    # if number of pixels available inferior number initial sample size
    if (Subset$nbPix2Sample < nb_Pix_To_Sample) {
      nb_Pix_To_Sample <- Subset$nbPix2Sample
      nb_partitions <- ceiling(nb_Pix_To_Sample / Pix_Per_Partition)
      Pix_Per_Partition <- floor(nb_Pix_To_Sample / nb_partitions)
      nb_Pix_To_Sample <- nb_partitions * Pix_Per_Partition
    }
    DataSubset <- Subset$DataSubset
    # # # assume that 1st data cleaning is enough...
    ## Uncommented June 5, 2019
    # clean reflectance data from inf and constant values
    CleanData <- rm_invariant_bands(DataSubset, SpectralFilter)
    DataSubset <- CleanData$DataMatrix
    SpectralFilter <- CleanData$Spectral
    print("perform PCA#2 on the subset image")
    if (TypePCA == "PCA" | TypePCA == "SPCA") {
      PCA_model <- pca(DataSubset, TypePCA)
    # } else if (TypePCA == "NLPCA") {
    #   print("performing NL-PCA with autoencoder")
    #   tic()
    #   PCA_model <- nlpca(DataSubset)
    #   toc()
    }
  }
  # Number of PCs computed and written in the PCA file: 30 if hyperspectral
  Nb_PCs <- dim(PCA_model$x)[2]
  if (Nb_PCs > NbPCs_To_Keep){
    Nb_PCs <- NbPCs_To_Keep
  }
  PCA_model$Nb_PCs <- Nb_PCs
  PCA_model$x <- NULL
  # CREATE PCA FILE CONTAINING ONLY SELECTED PCs
  print("Apply PCA model to the whole image")
  Output_Dir_PCA <- define_output_subdir(Output_Dir, Input_Image_File, TypePCA, "PCA")
  PCA_Files <- paste(Output_Dir_PCA, "OutputPCA_", Nb_PCs, "_PCs", sep = "")
  write_PCA_raster(Input_Image_File = Input_Image_File, Input_Mask_File = Input_Mask_File,
                   PCA_Path = PCA_Files, PCA_model = PCA_model, Spectral = SpectralFilter,
                   Nb_PCs = Nb_PCs, Continuum_Removal = Continuum_Removal, TypePCA = TypePCA,
                   nbCPU = nbCPU, MaxRAM = MaxRAM)
  # save workspace for this stage
  WS_Save <- paste(Output_Dir_PCA, "PCA_Info.RData", sep = "")
  my_list <- list("PCA_Files" = PCA_Files,"Pix_Per_Partition" =Pix_Per_Partition, "nb_partitions" = nb_partitions,
                  "MaskPath" = Input_Mask_File, "PCA_model" =PCA_model,"SpectralFilter"=   SpectralFilter)
  MaskPath = Input_Mask_File
  save(PCA_Files,Pix_Per_Partition, nb_partitions, MaskPath, PCA_model, SpectralFilter,file = WS_Save)
  return(my_list)
}

# perform filtering based on extreme values PCA identified through PCA
#
# @param Input_Image_File character. Path of the image to be processed
# @param HDR character. Path of the header file corresponding to the image to be processed
# @param Input_Mask_File character. Path of the mask raster corresponding to the image (keeps pixels = 1)
# @param Shade_Update character. Path of the updated mask raster corresponding to the image (keeps pixels = 1)
# @param Spectral list. spectral information from data
# @param Continuum_Removal boolean. set TRUE if continuum removal should be applied
# @param PCA_model dataframe. general parameters of the PCA
# @param PCsel numeric. PCs used to filter out extreme values
# @param TypePCA character. Set to PCA, SPCA or MNF
# @param nbCPU numeric. number of CPUs to be used in parallel
# @param MaxRAM numeric. indicator of RAM to be used to read image file
#
# @return Shade_Update = updated shade mask
# ' @importFrom stats sd
# ' @importFrom matlab ones
filter_PCA <- function(Input_Image_File, HDR, Input_Mask_File, Shade_Update,
                       Spectral, Continuum_Removal, PCA_model, PCsel, TypePCA,
                       nbCPU = 1, MaxRAM = 0.25) {

  # 1- get extreme values falling outside of mean +- 3SD for PCsel first components
  # compute mean and sd of the 5 first components of the sampled data
  MeanSub <- colMeans(PCA_model$x)
  SDSub <- apply(PCA_model$x, 2, sd)
  MinPC <- MeanSub - 3.0 * SDSub
  MaxPC <- MeanSub + 3.0 * SDSub

  # 2- update shade mask based on PCA values
  # 2.1- create hdr and binary files corresponding to updated mask
  HDR_Shade <- HDR
  HDR_Shade$description <- "Mask produced from PCA outlier filtering"
  HDR_Shade$bands <- 1
  HDR_Shade$`data type` <- 1
  HDR_Shade$`band names` <- "{Mask_PCA}"
  HDR_Shade$wavelength <- NULL
  HDR_Shade$fwhm <- NULL
  HDR_Shade$resolution <- NULL
  HDR_Shade$bandwidth <- NULL
  HDR_Shade$purpose <- NULL
  HDR_Shade$`default stretch` <- '0 1 linear'
  HDR_Shade$`default bands` <- NULL
  HDR_Shade$`data gain values` <- NULL
  HDR_Shade$`byte order` <- get_byte_order()
  headerFpath <- paste(Shade_Update, ".hdr", sep = "")
  write_ENVI_header(HDR_Shade, headerFpath)
  # create updated shade mask
  fidShade_Update <- file(
    description = Shade_Update, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidShade_Update)

  # 2.2- read image file sequentially
  Image_Format <- ENVI_type2bytes(HDR)
  Shade_Format <- ENVI_type2bytes(HDR_Shade)
  lenTot <- as.double(HDR$samples) * as.double(HDR$lines) * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image_Format$Bytes) / (1024^3)
  # maximum image size read at once. If image larger, then reads in multiple pieces
  LimitSizeGb <- MaxRAM
  if (ImSizeGb < LimitSizeGb) {
    Lines_Per_Read <- HDR$lines
    nbPieces <- 1
  } else {
    # nb of lines corresponding to LimitSizeGb
    OneLine <- as.double(HDR$samples) * as.double(HDR$bands) * Image_Format$Bytes
    Lines_Per_Read <- floor(LimitSizeGb * (1024^3) / OneLine)
    # number of pieces to split the image into
    nbPieces <- ceiling(HDR$lines / Lines_Per_Read)
  }
  # prepare for sequential processing: SeqRead_Image informs about byte location to read
  SeqRead_Image <- where_to_read(HDR, nbPieces)
  HDR_Shade <- read_ENVI_header(headerFpath)
  SeqRead_Shade <- where_to_read(HDR_Shade, nbPieces)
  Image_Format <- ENVI_type2bytes(HDR)

  print("Perform PCA on image subsets and filter data")
  # for each piece of image
  for (i in 1:nbPieces) {
    print(paste("PCA Piece #", i, "/", nbPieces))
    # read image and mask data
    nbLines <- SeqRead_Image$Lines_Per_Chunk[i]
    ImgFormat <- "2D"
    Image_Chunk <- read_image_subset(ImPath = Input_Image_File, HDR = HDR,
                                     Line_Start = SeqRead_Image$Line_Start[i],Lines_To_Read = nbLines,
                                     ImgFormat = ImgFormat)
    ImgFormat <- "Shade"
    if ((!Input_Mask_File == FALSE) & (!Input_Mask_File == "")) {
      Shade_Chunk <- read_image_subset(ImPath = Input_Mask_File, HDR = HDR_Shade,
                                       Line_Start = SeqRead_Image$Line_Start[i],Lines_To_Read = nbLines,
                                       ImgFormat = ImgFormat)
    } else {
      Shade_Chunk <- ones(nbLines * HDR$samples, 1)
    }
    keepShade <- which(Shade_Chunk == 1)
    Image_Chunk <- Image_Chunk[keepShade,]

    # apply Continuum removal if needed
    if (Continuum_Removal) {
      Image_Chunk <- apply_continuum_removal(Image_Chunk, Spectral, nbCPU = nbCPU)
    } else {
      if (!length(Spectral$WaterVapor) == 0) {
        Image_Chunk <- Image_Chunk[, -Spectral$WaterVapor]
      }
    }
    # remove constant bands if needed
    if (!length(Spectral$BandsNoVar) == 0) {
      Image_Chunk <- Image_Chunk[, -Spectral$BandsNoVar]
    }
    # Apply PCA
    if (TypePCA == "PCA" | TypePCA == "SPCA" | TypePCA == "MNF") {
      Image_Chunk <- scale(Image_Chunk, PCA_model$center, PCA_model$scale) %*% PCA_model$rotation[, PCsel]
    }

    # get PCA of the group of line and rearrange the data to write it correctly in the output file
    linetmp <- matrix(NA, ncol = ncol(Image_Chunk), nrow = (HDR$samples * nbLines))
    if (length(keepShade) > 0) {
      linetmp[keepShade, ] <- Image_Chunk
    }
    # find pixels which show extreme PC values
    ElimList <- list()
    for (pc in PCsel) {
      el0 <- matrix(which(linetmp[, pc] < MinPC[pc] | linetmp[, pc] > MaxPC[pc]), ncol = 1)
      if (length(el0) > 0) {
        ElimList <- c(ElimList, el0)
      }
    }
    elim <- unique(do.call("rbind", ElimList))
    if (length(elim) > 0) {
      Shade_Chunk[elim] <- 0
    }
    # files to write in
    fidOUT <- file(
      description = Shade_Update, open = "r+b", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    if (!SeqRead_Shade$ReadByte_Start[i] == 1) {
      seek(fidOUT, where = SeqRead_Shade$ReadByte_Start[i] - 1, origin = "start", rw = "write")
    }
    Shade_Chunk <- array(Shade_Chunk, c(nbLines, HDR_Shade$samples, 1))
    Shade_Chunk <- aperm(Shade_Chunk, c(2, 3, 1))
    writeBin(c(as.integer(Shade_Chunk)), fidOUT, size = 1, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
  }
  gc()
  return(Shade_Update)
}

# writes an ENVI image corresponding to PCA
#
# @param Input_Image_File path for the raster on which PCA is applied
# @param Input_Mask_File path for the corresponding mask
# @param PCA_Path path for resulting PCA
# @param PCA_model PCA model description
# @param Spectral spectral information to be used in the image
# @param Nb_PCs number of components kept in the resulting PCA raster
# @param CR boolean. If TRUE continuum removal is performed.
# @param TypePCA PCA, SPCA, NLPCA
# @param nbCPU number of CPUs to process data
# @param MaxRAM max RAM when initial image is read (in Gb)
#
# @return None
write_PCA_raster <- function(Input_Image_File, Input_Mask_File, PCA_Path, PCA_model,
                             Spectral, Nb_PCs, Continuum_Removal, TypePCA, nbCPU = 1, MaxRAM = 0.25) {
  ImPathHDR <- get_HDR_name(Input_Image_File)
  HDR <- read_ENVI_header(ImPathHDR)
  if (is.character(Input_Mask_File) && (Input_Mask_File != "")) {
    ShadeHDR <- get_HDR_name(Input_Mask_File)
    HDR_Shade <- read_ENVI_header(ShadeHDR)
  } else {
    HDR_Shade <- FALSE
  }
  # 1- create hdr and binary files corresponding to PCA file
  HDR_PCA <- HDR
  HDR_PCA$bands <- Nb_PCs
  HDR_PCA$`data type` <- 4
  HDR_PCA$interleave <- 'BIL'
  HDR_PCA$`default bands` <- NULL
  HDR_PCA$`wavelength units` <- NULL
  HDR_PCA$`z plot titles` <- NULL
  HDR_PCA$`data gain values` <- NULL
  HDR_PCA$`band names` <- paste('PC', 1:Nb_PCs, collapse = ", ")
  HDR_PCA$wavelength <- NULL
  HDR_PCA$fwhm <- NULL
  HDR_PCA$resolution <- NULL
  HDR_PCA$bandwidth <- NULL
  HDR_PCA$purpose <- NULL
  HDR_PCA$`default stretch` <- NULL
  HDR_PCA$`byte order` <- get_byte_order()
  headerFpath <- paste(PCA_Path, ".hdr", sep = "")
  write_ENVI_header(HDR_PCA, headerFpath)
  # create updated shade mask
  fidPCA <- file(
    description = PCA_Path, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidPCA)
  # apply PCA to the image
  # read image file sequentially
  Image_Format <- ENVI_type2bytes(HDR)
  PCA_Format <- ENVI_type2bytes(HDR_PCA)

  if (typeof(HDR_Shade) == 'list') {
    Shade_Format <- ENVI_type2bytes(HDR_Shade)
  } else if (typeof(HDR_Shade) == 'logical'){
    Shade_Format <- FALSE
  }
  # if (!HDR_Shade == FALSE) {
  #   Shade_Format <- ENVI_type2bytes(HDR_Shade)
  # } else {
  #   Shade_Format <- FALSE
  # }
  # prepare for sequential processing: SeqRead_Image informs about byte location to read
  nbPieces <- split_image(HDR, LimitSizeGb = MaxRAM)
  SeqRead_Image <- where_to_read(HDR, nbPieces)
  SeqRead_Shade <- where_to_read(HDR_Shade, nbPieces)
  SeqRead_PCA <- where_to_read(HDR_PCA, nbPieces)

  # for each piece of image
  for (i in 1:nbPieces) {
    print(paste("PCA Piece #", i, "/", nbPieces))
    # read image and mask data
    nbLines <- SeqRead_Image$Lines_Per_Chunk[i]
    ImgFormat <- "2D"
    Image_Chunk <- read_image_subset(ImPath = Input_Image_File, HDR = HDR,
                                     Line_Start = SeqRead_Image$Line_Start[i],Lines_To_Read = nbLines,
                                     ImgFormat = ImgFormat)
    if (typeof(HDR_Shade) == 'list') {
      ImgFormat <- "Shade"
      Shade_Chunk <- read_image_subset(ImPath = Input_Mask_File, HDR = HDR_Shade,
                                       Line_Start = SeqRead_Image$Line_Start[i],Lines_To_Read = nbLines,
                                       ImgFormat = ImgFormat)
      keepShade <- which(Shade_Chunk == 1)
      Image_Chunk <- Image_Chunk[keepShade, ]
    } else {
      keepShade <- matrix(1,ncol = 1,nrow = nrow(Image_Chunk))
    }
    # apply Continuum removal if needed
    if (Continuum_Removal) {
      Image_Chunk <- apply_continuum_removal(Image_Chunk, Spectral, nbCPU = nbCPU)
      ## added June 5, 2019
      if (length(Spectral$BandsNoVar) > 0) {
        Image_Chunk <- Image_Chunk[, -Spectral$BandsNoVar]
      }
    } else {
      # Eliminate water vapor
      Image_Chunk <- Image_Chunk[, Spectral$Bands2Keep]
      ## added June 5, 2019
      if (length(Spectral$BandsNoVar) > 0) {
        Image_Chunk <- Image_Chunk[, -Spectral$BandsNoVar]
      }
    }

    # Apply PCA
    if (TypePCA == "PCA" | TypePCA == "SPCA" | TypePCA == "MNF") {
      Image_Chunk <- scale(Image_Chunk, PCA_model$center, PCA_model$scale) %*% PCA_model$rotation[, 1:Nb_PCs]
    }

    # get PCA of the group of line and rearrange the data to write it correctly in the output file
    PCA_Chunk <- matrix(NA, ncol = Nb_PCs, nrow = (HDR$samples * nbLines))
    if (length(keepShade) > 0) {
      PCA_Chunk[keepShade, ] <- Image_Chunk
    }
    # files to write in
    fidOUT <- file(
      description = PCA_Path, open = "r+b", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    if (!SeqRead_PCA$ReadByte_Start[i] == 1) {
      nbSkip <- (SeqRead_PCA$ReadByte_Start[i] - 1) * PCA_Format$Bytes
      seek(fidOUT, where = nbSkip, origin = "start", rw = "write")
    }
    PCA_Chunk <- array(PCA_Chunk, c(nbLines, HDR_PCA$samples, HDR_PCA$bands))
    PCA_Chunk <- aperm(PCA_Chunk, c(2, 3, 1))
    writeBin(c(PCA_Chunk), fidOUT, size = PCA_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
  }
  list <- ls()
  rm(list = list)
  gc()
  return(invisible())
}

# Function to perform PCA on a matrix
#
# @param X matrix to apply PCA on
# @param type PCA (no rescale) or SPCA (rescale)
#
# @return list of PCA parameters (PCs from X, mean, eigenvectors and values)
#' @importFrom stats prcomp
pca <- function(X, type) {
  p <- ncol(X)
  if (type == "SPCA") {
    modPCA <- prcomp(X, scale = TRUE)

  } else if (type == "PCA") {
    modPCA <- prcomp(X, scale = FALSE)

  }

  return(modPCA)
}

# defines the number of pixels per iteration based on a trade-off between image size and sample size per iteration
#
# @param ImNames Path and name of the images to be processed
# @param nb_partitions number of iterations peformed to average diversity indices
#
# @return Pix_Per_Partition number of pixels per iteration
define_pixels_per_iter <- function(ImNames, nb_partitions = nb_partitions) {
  Input_Image_File <- ImNames$Input.Image
  Input_Mask_File <- ImNames$Mask_list
  # define dimensions of the image
  ImPathHDR <- get_HDR_name(Input_Image_File)
  HDR <- read_ENVI_header(ImPathHDR)
  Image_Format <- ENVI_type2bytes(HDR)
  ipix <- as.double(HDR$lines)
  jpix <- as.double(HDR$samples)
  nbPixels <- ipix * jpix
  lenTot <- nbPixels * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image_Format$Bytes) / (1024^3)
  # if shade mask, update number of pixels
  if (is.character(Input_Mask_File) && (Input_Mask_File != "")) {
    # read shade mask
    fid <- file(
      description = Input_Mask_File, open = "rb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    ShadeMask <- readBin(fid, integer(), n = nbPixels, size = 1)
    close(fid)
    ShadeMask <- aperm(array(ShadeMask, dim = c(jpix, ipix)))
    # number of sunlit pixels
    nbPixels_Sunlit <- length(which(ShadeMask == 1))
  } else {
    nbPixels_Sunlit <- nbPixels
  }

  # adjust the number of pixels per iteration
  # trade-off between number of pixels and total pixel size
  # maximum number of pixels to be used
  Max_Pixel_Per_Iter <- 5000000/nb_partitions
  Max_Pixel_Per_Iter_Size <- nb_partitions * (Max_Pixel_Per_Iter * as.double(HDR$bands) * Image_Format$Bytes) / (1024^3)
  # maximum amount of data (Gb) to be used
  Max_Size_Per_Iter <- 0.3
  # amount of data available after first mask
  ImDataGb <- (nbPixels_Sunlit * as.double(HDR$bands) * Image_Format$Bytes) / (1024^3)

  # if Max_Pixel_Per_Iter correspond to reasonable size (<Max_Size_Per_Iter)
  if (Max_Pixel_Per_Iter_Size < Max_Size_Per_Iter) {
    # if enough data
    if (ImDataGb >= Max_Pixel_Per_Iter_Size) {
      Pix_Per_Partition <- Max_Pixel_Per_Iter
    } else if (ImDataGb < Max_Pixel_Per_Iter_Size) {
      # define number of pixels corresponding to ImDataGb
      Pix_Per_Partition <- floor(nbPixels_Sunlit / nb_partitions)
    }
    # if size too important, adjust number of pixels to match Max_Size_Per_Iter
  } else if (Max_Pixel_Per_Iter_Size >= Max_Size_Per_Iter) {
    Pix_Per_Partition <- floor((Max_Size_Per_Iter * (1024^3) / (as.double(HDR$bands) * Image_Format$Bytes)) / nb_partitions)
  }
  return(Pix_Per_Partition)
}

#' Check if principal components are properly selected as expected by the method
#'
#' @param Input_Image_File character. Path of the image to be processed
#' @param Output_Dir character. Path for output directory
#' @param PCA_Files character. Path of the PCA image
#' @param TypePCA character. Type of PCA: choose either "PCA" or "SPCA"
#' @param File_Open Boolean. Set to TRUE for file to open automatically
#'
#' @return Sel_PC
#' @importFrom utils file.edit
#' @export
select_PCA_components <- function(Input_Image_File, Output_Dir, PCA_Files, TypePCA = "SPCA", File_Open = FALSE) {
  message("")
  message("*********************************************************")
  message("Please check following PCA file:")
  print(PCA_Files)
  message("*********************************************************")
  Image_Name <- strsplit(basename(Input_Image_File), "\\.")[[1]][1]
  Output_Dir_Full <- paste(Output_Dir, "/", Image_Name, "/", TypePCA, "/", sep = "")
  Sel_PC <- paste(Output_Dir_Full, "/PCA/Selected_Components.txt", sep = "")
  message("list the principal components that will be used to estimate biodiversity in the file")
  message("")
  print(Sel_PC)
  message("")
  message("Then press ENTER")

  if (!file.exists(Sel_PC)) {
    file.create(Sel_PC)
  }
  if (File_Open == TRUE) {
    file.edit(Sel_PC, title=basename(Sel_PC),editor = "internal")
  }
  message("*********************************************************")
  message("")
  readline(prompt = "")
  return(Sel_PC)
}

#' this function performs rescaling and
#' either defines min and max from each feature in a data set,
#' or applies the transformation based on a previously defined min and max
#'
#' @param x numeric. data matrix
#' @param mode character. 'define' or 'apply'
#' @param MinX numeric. if 'apply'
#' @param MaxX numeric. if 'apply'
#'
#' @return rescaled data, min and max values
#' @export
minmax <- function(x, mode = "define", MinX = FALSE, MaxX = FALSE) {
  ## Function to
  if (mode == "define") {
    MinX <- min(x)
    MaxX <- max(x)
    (x - MinX) / (MaxX - MinX)
  } else if (mode == "apply") {
    (x - MinX) / (MaxX - MinX)
  }
  my_list <- list("data" = x, "MinX" = MinX, "MaxX" = MaxX)
  return(my_list)
}


#' this function
#'
#' @param X numeric. data matrix
#' @param coordPix numeric.
#'
#' @return rescaled data, min and max values
#' @importFrom matlab ndims
#' @export
# If coordPix is not NULL, X and coordPix are exepected to have the same order,
# i.e. coordPix[1, ] corresponds to X[1, ], coordPix[2, ] corresponds to X[2, ], ...
noise <- function(X, coordPix=NULL){
  if(is.null(coordPix)){
    if(ndims(X)!=3)
      stop('X is expected to be a 3D array: y,x,band for row,col,depth.')
    Xdim = dim(X)
    # Shift x/y difference
    Y = ((X[2:Xdim[1],2:Xdim[2],]-X[1:(Xdim[1]-1),2:Xdim[2],]) +
               (X[2:Xdim[1],2:Xdim[2],]-X[2:Xdim[1],1:(Xdim[2]-1),]))/2
  }else{
    if(!all(c('Kind', 'id') %in% colnames(coordPix)))
      stop("Columns 'Kind' and 'id' are missing in coordPix.")
    kernel = matrix(0, 3, 3)
    kernel[c(5, 6, 8)]=c(1, -1/2, -1/2)

    if(!identical(order(coordPix$id), 1:nrow(coordPix)))
      stop("coordPix is not ordered along column 'id'. Order coordPix as well as X before trying again.")

    Y=0
    for(ik in which(kernel!=0)){
      Y = Y + X[coordPix$Kind==ik,]*kernel[ik]
    }
  }
  return(Y)
}

#' XXX
#'
#' @param X numeric.
#' @param kernel numeric.
#'
#' @return coordinates of the pixels
#' @export
# used in noise
coordPix_kernel <- function(X, kernel){
  mesh = matlab::meshgrid(1:nrow(X), 1:ncol(X))
  Row <- as.numeric(mesh$y) # row
  Column <- as.numeric(mesh$x) # column

  coordPixK = list()
  mesh=matlab::meshgrid(-(ncol(kernel)%/%2):(ncol(kernel)%/%2), -(nrow(kernel)%/%2):(nrow(kernel)%/%2))
  for(p in which(kernel!=0)){
    coordPixK[[p]] = data.table(row = Row+mesh$y[p], col = Column+mesh$x[p], id=1:length(Row))
  }
  coordPix = rbindlist(coordPixK, idcol='Kind')
  # Order along coordPix$id for further use in noise, mnf
  setorder(coordPix, 'id')
  return(coordPix)
}

#' Function to perform MNF
#'
#' @param X numeric. matrix to apply MNF on
#' @param coordPix dataframe to compute noise, cf get_random_subset_from_image
#' @param retx boolean.
#'
#' @return results of MNF applied on matrix
#' @importFrom stats cov
#' @export
# used in noise

# TODO: faire 2 fonctions: mnf et mnf.subset, de même pour noise, noise.subset
mnf <- function(X, coordPix=NULL, retx=TRUE){
  if(any(is.na(X))){
    stop('Pixels with NAs found in X. Remove NA pixels before trying again.')
  }

  if(length(dim(X))>3)
    stop('X has more than 3 dimensions.')

  nz <- noise(X, coordPix)
  Xdim = dim(X)

  if(is.null(coordPix) && length(dim(X))>2){
    X = matrix(X[1:(Xdim[1]-1), 1:(Xdim[2]-1),], nrow = Xdim[1]*Xdim[2])
    nz = matrix(nz, nrow = Xdim[1]*Xdim[2])
  }

  Xc = scale(X, center = T, scale = F)

  covNoise = cov(nz)
  covXc = cov(Xc)
  eig = eigen(solve(covNoise)%*%covXc)
  colnames(eig$vectors) = paste0('PC', 1:ncol(eig$vectors))
  modMNF = list(sdev = sqrt(eig$values), rotation = eig$vectors,
                center = colMeans(X), scale = FALSE)
  attr(modMNF, 'class') <- 'prcomp'
  #   eig_pairs = tofsims:::EigenDecompose(covXc, covNoise, 1, nrow(covNoise))
  #   vord = order(Re(eig_pairs$eigval), decreasing = T)
  #   eig_pairs$eigval = Re(eig_pairs$eigval)[vord]
  #   eig_pairs$eigvec = Re(eig_pairs$eigvec[, vord])
  #   modMNF = list(rotation=eig_pairs$eigvec,
  #                 sdev=sqrt(eig_pairs$eigval),
  #                 center=colMeans(X),
  #                 scale=FALSE)
  if(retx==T)
    modMNF$x= array(Xc %*% modMNF$rotation, dim = Xdim)

  return(modMNF)
}
