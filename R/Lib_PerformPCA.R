# ==============================================================================
# biodivMapR
# Lib_PerformPCA.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
# ==============================================================================
# This Library is used to perform PCA on raster prior to diversity mapping
# ==============================================================================

#' Performs PCA for all images and create PCA file with either all or a selection of PCs
#'
#' @param ImNames Path and name of the images to be processed
#' @param Output.Dir output directory
#' @param Continuum.Removal boolean: should continuum removal be applied?
#' @param TypePCA Type of PCA (PCA, SPCA, NLPCA...)
#' @param FilterPCA boolean. If TRUE 2nd filtering based on PCA
#' @param Excluded.WL boolean. Water Vapor Absorption domains (in nanometers). Can also be used to exclude spectific domains
#' @param NbIter numeric. Number of iteration to estimate diversity from the raster.
#' @param nbCPU numeric. Number fo CPUs in use.
#' @param MaxRAM numeric. Maximum size of chunk in GB to limit RAM allocation when reading image file.
#'
#' @return list of paths corresponding to resulting PCA files
#' @export
Perform.PCA.Image <- function(ImPath, ImPathShade, Output.Dir, Continuum.Removal = TRUE, TypePCA = "SPCA", FilterPCA = FALSE, Excluded.WL = FALSE, NbIter = 20, nbCPU = 1, MaxRAM = 0.25) {
  # define the path corresponding to image, mask and output directory
  ImNames <- list()
  ImNames$Input.Image <- ImPath
  ImNames$Mask.list <- ImPathShade
  Output.Dir2 <- Define.Output.Dir(Output.Dir, ImPath, TypePCA)
  # Identify water vapor absorption bands in image and possibly other spectral domains to discard
  Spectral <- Exclude.Spectral.Domains(ImPath, Excluded.WL = Excluded.WL)
  # Extract valid data subset and check validity
  print("Extract pixels from the images to perform PCA on a subset")
  # define number of pixels to be extracted from the image for each iteration
  Pix.Per.Iter <- Define.Pixels.Per.Iter(ImNames, NbIter = NbIter)
  nb.Pix.To.Sample <- NbIter * Pix.Per.Iter
  ImPathHDR <- Get.HDR.Name(ImPath)
  HDR <- read.ENVI.header(ImPathHDR)
  # extract a random selection of pixels from image
  Subset <- Get.Random.Subset.From.Image(ImPath, HDR, ImPathShade, NbIter, Pix.Per.Iter)
  # if needed, apply continuum removal
  if (Continuum.Removal == TRUE) {
    Subset$DataSubset <- Apply.Continuum.Removal(Subset$DataSubset, Spectral, nbCPU = nbCPU)
  }
  # if number of pixels available inferior number initial sample size
  if (Subset$nbPix2Sample < nb.Pix.To.Sample) {
    nb.Pix.To.Sample <- Subset$nbPix2Sample
    NbIter <- ceiling(nb.Pix.To.Sample / Pix.Per.Iter)
    Pix.Per.Iter <- floor(nb.Pix.To.Sample / NbIter)
    nb.Pix.To.Sample <- NbIter * Pix.Per.Iter
  }
  DataSubset <- Subset$DataSubset
  # clean reflectance data from inf and constant values
  CleanData <- Check.Data(DataSubset, Spectral)
  DataSubset <- CleanData$DataMatrix
  Spectral <- CleanData$Spectral

  # Compute PCA #1 on DataSubset
  print("perform PCA#1 on the subset image")
  if (TypePCA == "PCA" | TypePCA == "SPCA") {
    PCA.model <- pca(DataSubset, TypePCA)
  } else if (TypePCA == "NLPCA") {
    print("performing NL-PCA with autoencoder")
    print("Make sure you properly installed and defined python environment if using this functionality")
    tic()
    PCA.model <- nlpca(DataSubset)
    toc()
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
    if (dim(PCA.model$dataPCA)[2] > 5) {
      PCsel <- 1:5
    } else {
      PCsel <- 1:dim(PCA.model$dataPCA)[2]
    }
    Shade.Update <- paste(Output.Dir2, "ShadeMask_Update_PCA", sep = "")
    Subset.PCA <- PCA.model$dataPCA[, PCsel]
    ImPathShade <- Filter.PCA(ImPath, HDR, ImPathShade, Shade.Update, Spectral, Continuum.Removal, PCA.model, PCsel, TypePCA, nbCPU = nbCPU, MaxRAM = MaxRAM)
    ## Compute PCA 2 based on the updated shade mask ##
    # extract a random selection of pixels from image
    Subset <- Get.Random.Subset.From.Image(ImPath, HDR, ImPathShade, NbIter, Pix.Per.Iter)
    # if needed, apply continuum removal
    if (Continuum.Removal == TRUE) {
      Subset$DataSubset <- Apply.Continuum.Removal(Subset$DataSubset, Spectral, nbCPU = nbCPU)
    }
    # if number of pixels available inferior number initial sample size
    if (Subset$nbPix2Sample < nb.Pix.To.Sample) {
      nb.Pix.To.Sample <- Subset$nbPix2Sample
      NbIter <- ceiling(nb.Pix.To.Sample / Pix.Per.Iter)
      Pix.Per.Iter <- floor(nb.Pix.To.Sample / NbIter)
      nb.Pix.To.Sample <- NbIter * Pix.Per.Iter
    }
    DataSubset <- Subset$DataSubset
    # # # assume that 1st data cleaning is enough...
    ## Uncommented June 5, 2019
    # clean reflectance data from inf and constant values
    CleanData <- Check.Data(DataSubset, Spectral)
    DataSubset <- CleanData$DataMatrix
    Spectral <- CleanData$Spectral
    print("perform PCA#2 on the subset image")
    if (TypePCA == "PCA" | TypePCA == "SPCA") {
      PCA.model <- pca(DataSubset, TypePCA)
    } else if (TypePCA == "NLPCA") {
      print("performing NL-PCA with autoencoder")
      tic()
      PCA.model <- nlpca(DataSubset)
      toc()
    }
  }
  # Number of PCs computed and written in the PCA file: 30 if hyperspectral
  Nb.PCs <- dim(PCA.model$dataPCA)[2]
  if (Nb.PCs > 30) {
    Nb.PCs <- 30
  }
  # CREATE PCA FILE CONTAINING ONLY SELECTED PCs
  print("Apply PCA model to the whole image")
  Output.Dir.PCA <- Define.Output.SubDir(Output.Dir, ImPath, TypePCA, "PCA")
  PCA.Path <- paste(Output.Dir.PCA, "OutputPCA_", Nb.PCs, "_PCs", sep = "")
  Create.PCA.Image(ImPath, ImPathShade, PCA.Path, PCA.model, Spectral, Nb.PCs, Continuum.Removal, TypePCA, nbCPU, MaxRAM = MaxRAM)
  # save workspace for this stage
  WS_Save <- paste(Output.Dir2, "PCA_Info.RData", sep = "")
  save(PCA.model, Pix.Per.Iter, NbIter, ImPathShade, file = WS_Save)
  PCA.Files <- PCA.Path
  return(PCA.Files)
}

# perform filtering based on extreme values PCA identified through PCA
#
# @param ImPath image path
# @param HDR hdr file file correspmonding to ImPath
# @param ImPathShade shade file path
# @param Shade.Update updated shade mask
# @param Spectral spectral information from data
# @param CR logical: does continuum removal need to be performed?
# @param PCA.model general parameters of the PCA
# @param TypePCA
# @param nbCPU
# @param MaxRAM
# @param PCsel PCs used to filter out extreme values
#
# @return Shade.Update = updated shade mask
Filter.PCA <- function(ImPath, HDR, ImPathShade, Shade.Update, Spectral, CR, PCA.model, PCsel, TypePCA, nbCPU = 1, MaxRAM = 0.25) {

  # 1- get extreme values falling outside of mean +- 3SD for PCsel first components
  # compute mean and sd of the 5 first components of the sampled data
  MeanSub <- colMeans(PCA.model$dataPCA)
  SDSub <- apply(PCA.model$dataPCA, 2, sd)
  MinPC <- MeanSub - 3.0 * SDSub
  MaxPC <- MeanSub + 3.0 * SDSub

  # 2- update shade mask based on PCA values
  # 2.1- create hdr and binary files corresponding to updated mask
  HDR.Shade <- HDR
  HDR.Shade$bands <- 1
  HDR.Shade$`data type` <- 1
  HDR.Shade$`band names` <- "{Mask_PCA}"
  HDR.Shade$wavelength <- NULL
  HDR.Shade$fwhm <- NULL
  HDR.Shade$resolution <- NULL
  HDR.Shade$bandwidth <- NULL
  HDR.Shade$purpose <- NULL
  HDR.Shade$`byte order` <- Get.Byte.Order()
  headerFpath <- paste(Shade.Update, ".hdr", sep = "")
  write.ENVI.header(HDR.Shade, headerFpath)
  # create updated shade mask
  fidShade.Update <- file(
    description = Shade.Update, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidShade.Update)

  # 2.2- read image file sequentially
  Image.Format <- ENVI.Type2Bytes(HDR)
  Shade.Format <- ENVI.Type2Bytes(HDR.Shade)
  lenTot <- as.double(HDR$samples) * as.double(HDR$lines) * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image.Format$Bytes) / (1024^3)
  # maximum image size read at once. If image larger, then reads in multiple pieces
  LimitSizeGb <- MaxRAM
  if (ImSizeGb < LimitSizeGb) {
    Lines.Per.Read <- HDR$lines
    nbPieces <- 1
  } else {
    # nb of lines corresponding to LimitSizeGb
    OneLine <- as.double(HDR$samples) * as.double(HDR$bands) * Image.Format$Bytes
    Lines.Per.Read <- floor(LimitSizeGb * (1024^3) / OneLine)
    # number of pieces to split the image into
    nbPieces <- ceiling(HDR$lines / Lines.Per.Read)
  }
  # prepare for sequential processing: SeqRead.Image informs about byte location to read
  SeqRead.Image <- Where.To.Read(HDR, nbPieces)
  HDR.Shade <- read.ENVI.header(headerFpath)
  SeqRead.Shade <- Where.To.Read(HDR.Shade, nbPieces)
  Image.Format <- ENVI.Type2Bytes(HDR)

  print("Perform PCA on image subsets and filter data")
  # for each piece of image
  for (i in 1:nbPieces) {
    print(paste("PCA Piece #", i, "/", nbPieces))
    # read image and mask data
    Byte.Start <- SeqRead.Image$ReadByte.Start[i]
    nbLines <- SeqRead.Image$Lines.Per.Chunk[i]
    lenBin <- SeqRead.Image$ReadByte.End[i] - SeqRead.Image$ReadByte.Start[i] + 1
    ImgFormat <- "2D"
    Image.Chunk <- Read.Image.Subset(ImPath, HDR, Byte.Start, lenBin, nbLines, Image.Format, ImgFormat)

    Byte.Start <- SeqRead.Shade$ReadByte.Start[i]
    nbLines <- SeqRead.Shade$Lines.Per.Chunk[i]
    lenBin <- SeqRead.Shade$ReadByte.End[i] - SeqRead.Shade$ReadByte.Start[i] + 1
    ImgFormat <- "Shade"
    if (!ImPathShade == "") {
      Shade.Chunk <- Read.Image.Subset(ImPathShade, HDR.Shade, Byte.Start, lenBin, nbLines, Shade.Format, ImgFormat)
    } else {
      Shade.Chunk <- ones(nbLines * HDR$samples, 1)
    }

    elimShade <- which(Shade.Chunk == 0)
    keepShade <- which(Shade.Chunk == 1)
    Image.Chunk <- Image.Chunk[keepShade, ]

    # these two lines are perfomed in Apply.Continuum.Removal
    # # Eliminate water vapor
    # Image.Chunk     = Image.Chunk[,Spectral$Bands2Keep]
    # apply Continuum removal if needed
    if (CR == TRUE) {
      Image.Chunk <- Apply.Continuum.Removal(Image.Chunk, Spectral, nbCPU = nbCPU)
    }
    # remove constant bands if needed
    if (!length(Spectral$BandsNoVar) == 0) {
      Image.Chunk <- Image.Chunk[, -Spectral$BandsNoVar]
    }
    # Apply PCA
    if (TypePCA == "PCA" | TypePCA == "SPCA") {
      Image.Chunk <- t(t(PCA.model$eiV[, PCsel]) %*% t(Center.Reduce(Image.Chunk, PCA.model$mu, PCA.model$scale)))
    } else if (TypePCA == "NLPCA") {
      # apply a transformation to the dataset
      Image.Chunk <- apply(Image.Chunk, 2, minmax, mode = "apply", MinX = PCA.model$MinVal, MaxX = PCA.model$MaxVal)
      Image.Chunk2 <- matrix(NA, nrow = length(Image.Chunk[[1]]$data), ncol = length(Image.Chunk))
      for (ii in 1:length(Image.Chunk)) {
        Image.Chunk2[, ii] <- matrix(Image.Chunk[[ii]]$data, ncol = 1)
      }
      Image.Chunk <- Image.Chunk2
      rm(Image.Chunk2)
      intermediate_layer_model <- keras_model(inputs = PCA.model$Model$input, outputs = get_layer(PCA.model$Model, "bottleneck")$output)

      # # use multithread to apply autoencoder...
      # plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
      # Image.Chunk     = splitRows(Image.Chunk, nbCPU)
      # Image.Chunk     = future_lapply(Image.Chunk,FUN = predict.from.NLPCA,model  = intermediate_layer_model,future.scheduling = 1.0)
      # plan(sequential)
      # Image.Chunk     = rbind(Image.Chunk)
      Image.Chunk <- predict(intermediate_layer_model, Image.Chunk)
    }

    # get PCA of the group of line and rearrange the data to write it correctly in the output file
    linetmp <- matrix(NA, ncol = ncol(Image.Chunk), nrow = (HDR$samples * HDR$lines))
    if (length(keepShade) > 0) {
      linetmp[keepShade, ] <- Image.Chunk
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
      Shade.Chunk[elim] <- 0
    }
    # files to write in
    fidOUT <- file(
      description = Shade.Update, open = "r+b", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    if (!SeqRead.Shade$ReadByte.Start[i] == 1) {
      seek(fidOUT, where = SeqRead.Shade$ReadByte.Start[i] - 1, origin = "start", rw = "write")
    }
    Shade.Chunk <- array(Shade.Chunk, c(nbLines, HDR.Shade$samples, 1))
    Shade.Chunk <- aperm(Shade.Chunk, c(2, 3, 1))
    writeBin(c(as.integer(Shade.Chunk)), fidOUT, size = 1, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
  }
  return(Shade.Update)
}

# writes an ENVI image corresponding to PCA
#
# @param ImPath path for the raster on which PCA is applied
# @param ImPathShade path for the corresponding shade mask
# @param PCA.Path path for resulting PCA
# @param PCA.model PCA model description
# @param Spectral spectral information to be used in the image
# @param Nb.PCs number of components kept in the resulting PCA raster
# @param CR Shoudl continuum removal be performed?
# @param TypePCA PCA, SPCA, NLPCA
# @param nbCPU number of CPUs to process data
# @param MaxRAM max RAM when initial image is read (in Gb)
#
# @return
Create.PCA.Image <- function(ImPath, ImPathShade, PCA.Path, PCA.model, Spectral, Nb.PCs, CR, TypePCA, nbCPU = 1, MaxRAM = 0.25) {
  ImPathHDR <- Get.HDR.Name(ImPath)
  HDR <- read.ENVI.header(ImPathHDR)
  ShadeHDR <- Get.HDR.Name(ImPathShade)
  HDR.Shade <- read.ENVI.header(ShadeHDR)
  # 1- create hdr and binary files corresponding to PCA file
  HDR.PCA <- HDR
  HDR.PCA$bands <- Nb.PCs
  HDR.PCA$`data type` <- 4
  PCs <- list()
  for (i in 1:Nb.PCs) {
    PCs <- c(PCs, paste("PC", i))
  }
  PCs <- paste(PCs, collapse = ", ")
  HDR.PCA$`band names` <- PCs
  HDR.PCA$wavelength <- NULL
  HDR.PCA$fwhm <- NULL
  HDR.PCA$resolution <- NULL
  HDR.PCA$bandwidth <- NULL
  HDR.PCA$purpose <- NULL
  HDR.PCA$`default stretch` <- NULL
  HDR.PCA$`byte order` <- Get.Byte.Order()
  headerFpath <- paste(PCA.Path, ".hdr", sep = "")
  write.ENVI.header(HDR.PCA, headerFpath)
  # create updated shade mask
  fidPCA <- file(
    description = PCA.Path, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidPCA)
  # apply PCA to the image
  # read image file sequentially
  Image.Format <- ENVI.Type2Bytes(HDR)
  PCA.Format <- ENVI.Type2Bytes(HDR.PCA)
  Shade.Format <- ENVI.Type2Bytes(HDR.Shade)
  # prepare for sequential processing: SeqRead.Image informs about byte location to read
  nbPieces <- Split.Image(HDR, LimitSizeGb = MaxRAM)
  SeqRead.Image <- Where.To.Read(HDR, nbPieces)
  SeqRead.Shade <- Where.To.Read(HDR.Shade, nbPieces)
  SeqRead.PCA <- Where.To.Read(HDR.PCA, nbPieces)

  # for each piece of image
  for (i in 1:nbPieces) {
    print(paste("PCA Piece #", i, "/", nbPieces))
    # read image and mask data
    Byte.Start <- SeqRead.Image$ReadByte.Start[i]
    nbLines <- SeqRead.Image$Lines.Per.Chunk[i]
    lenBin <- SeqRead.Image$ReadByte.End[i] - SeqRead.Image$ReadByte.Start[i] + 1
    ImgFormat <- "2D"
    Image.Chunk <- Read.Image.Subset(ImPath, HDR, Byte.Start, lenBin, nbLines, Image.Format, ImgFormat)

    Byte.Start <- SeqRead.Shade$ReadByte.Start[i]
    nbLines <- SeqRead.Shade$Lines.Per.Chunk[i]
    lenBin <- SeqRead.Shade$ReadByte.End[i] - SeqRead.Shade$ReadByte.Start[i] + 1
    ImgFormat <- "Shade"
    Shade.Chunk <- Read.Image.Subset(ImPathShade, HDR.Shade, Byte.Start, lenBin, nbLines, Shade.Format, ImgFormat)

    elimShade <- which(Shade.Chunk == 0)
    keepShade <- which(Shade.Chunk == 1)
    Image.Chunk <- Image.Chunk[keepShade, ]

    # apply Continuum removal if needed
    if (CR == TRUE) {
      Image.Chunk <- Apply.Continuum.Removal(Image.Chunk, Spectral, nbCPU = nbCPU)
      ## added June 5, 2019
      if (length(Spectral$BandsNoVar) > 0) {
        Image.Chunk <- Image.Chunk[, -Spectral$BandsNoVar]
      }
    } else {
      # Eliminate water vapor
      Image.Chunk <- Image.Chunk[, Spectral$Bands2Keep]
      ## added June 5, 2019
      if (length(Spectral$BandsNoVar) > 0) {
        Image.Chunk <- Image.Chunk[, -Spectral$BandsNoVar]
      }
    }

    # !!!! EDIT 06-Feb-2019: comment these lines
    # # remove constant bands if needed
    # if (!length(Spectral$BandsNoVar)==0) {
    #   Image.Chunk   = Image.Chunk[,-Spectral$BandsNoVar]
    # }
    # Apply PCA
    # Apply PCA
    if (TypePCA == "PCA" | TypePCA == "SPCA") {
      Image.Chunk <- t(t(PCA.model$eiV[, 1:Nb.PCs]) %*% t(Center.Reduce(Image.Chunk, PCA.model$mu, PCA.model$scale)))
    } else if (TypePCA == "NLPCA") {
      # apply a transformation to the dataset
      Image.Chunk <- apply(Image.Chunk, 2, minmax, mode = "apply", MinX = PCA.model$MinVal, MaxX = PCA.model$MaxVal)
      Image.Chunk2 <- matrix(NA, nrow = length(Image.Chunk[[1]]$data), ncol = length(Image.Chunk))
      for (ii in 1:length(Image.Chunk)) {
        Image.Chunk2[, ii] <- matrix(Image.Chunk[[ii]]$data, ncol = 1)
      }
      Image.Chunk <- Image.Chunk2
      rm(Image.Chunk2)
      intermediate_layer_model <- keras_model(inputs = PCA.model$Model$input, outputs = get_layer(PCA.model$Model, "bottleneck")$output)

      # # use multithread to apply autoencoder...
      # plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
      # Image.Chunk     = splitRows(Image.Chunk, nbCPU)
      # Image.Chunk     = future_lapply(Image.Chunk,FUN = predict.from.NLPCA,model  = intermediate_layer_model,future.scheduling = 1.0)
      # plan(sequential)
      # Image.Chunk     = rbind(Image.Chunk)
      Image.Chunk <- predict(intermediate_layer_model, Image.Chunk)
    }

    # get PCA of the group of line and rearrange the data to write it correctly in the output file
    PCA.Chunk <- matrix(NA, ncol = Nb.PCs, nrow = (HDR$samples * nbLines))
    if (length(keepShade) > 0) {
      PCA.Chunk[keepShade, ] <- Image.Chunk
    }
    # files to write in
    fidOUT <- file(
      description = PCA.Path, open = "r+b", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    if (!SeqRead.PCA$ReadByte.Start[i] == 1) {
      # nbSkip      = (as.double(SeqRead.PCA$ReadByte.Start[i])*as.double(PCA.Format$Bytes))-1
      # nbSkip      = (as.double(SeqRead.PCA$ReadByte.Start[i])*as.double(PCA.Format$Bytes))-as.double(1)
      nbSkip <- (SeqRead.PCA$ReadByte.Start[i] - 1) * PCA.Format$Bytes
      seek(fidOUT, where = nbSkip, origin = "start", rw = "write")
    }
    PCA.Chunk <- array(PCA.Chunk, c(nbLines, HDR.PCA$samples, HDR.PCA$bands))
    PCA.Chunk <- aperm(PCA.Chunk, c(2, 3, 1))
    # writeBin(as.numeric(c(PCA.Chunk)), fidOUT, size = PCA.Format$Bytes,endian = .Platform$endian)
    writeBin(c(PCA.Chunk), fidOUT, size = PCA.Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    close(fidOUT)
  }
  list <- ls()
  rm(list = list)
  gc()
  return()
}


# Function to perform NLPCA
#
# @param DataSubset matrix to apply NLPCA on
#
# @return list of PCA parameters (PCs from X, mean, eigenvectors and values)
nlpca <- function(DataSubset) {

  # put between 0 and 1
  x_train <- apply(DataSubset, 2, minmax, mode = "define")
  Subset <- list()
  nb.Vars <- ncol(DataSubset)
  Subset$DataSubset <- matrix(NA, nrow = nrow(DataSubset), ncol = nb.Vars)
  for (i in 1:nb.Vars) {
    Subset$DataSubset[, i] <- matrix(x_train[[i]]$data, ncol = 1)
    Subset$MinVal[i] <- x_train[[i]]$MinX
    Subset$MaxVal[i] <- x_train[[i]]$MaxX
  }
  # define number of pixels to be used for NL-PCA
  nbSubsamples.NLPCA <- 100000
  # autoencoder in keras
  # set training data
  x_train <- as.matrix(Subset$DataSubset)

  # set model
  model <- keras_model_sequential()
  if (nb.Vars < 12) {
    model %>%
      layer_dense(units = 6, activation = "tanh", input_shape = ncol(x_train)) %>%
      layer_dense(units = 3, activation = "tanh", name = "bottleneck") %>%
      layer_dense(units = 6, activation = "tanh") %>%
      layer_dense(units = ncol(x_train))
  } else if (nb.Vars > 100) {
    model %>%
      layer_dense(units = 100, activation = "tanh", input_shape = ncol(x_train)) %>%
      layer_dense(units = 80, activation = "tanh") %>%
      layer_dense(units = 60, activation = "tanh") %>%
      layer_dense(units = 40, activation = "tanh") %>%
      layer_dense(units = 20, activation = "tanh", name = "bottleneck") %>%
      layer_dense(units = 40, activation = "tanh") %>%
      layer_dense(units = 60, activation = "tanh") %>%
      layer_dense(units = 80, activation = "tanh") %>%
      layer_dense(units = 100, activation = "tanh") %>%
      layer_dense(units = ncol(x_train))
  } else if (nb.Vars <= 100) {
    model %>%
      layer_dense(units = 80, activation = "tanh", input_shape = ncol(x_train)) %>%
      layer_dense(units = 60, activation = "tanh") %>%
      layer_dense(units = 40, activation = "tanh") %>%
      layer_dense(units = 20, activation = "tanh", name = "bottleneck") %>%
      layer_dense(units = 40, activation = "tanh") %>%
      layer_dense(units = 60, activation = "tanh") %>%
      layer_dense(units = 80, activation = "tanh") %>%
      layer_dense(units = ncol(x_train))
  }
  # view model layers
  summary(model)
  # compile model
  model %>% keras::compile(
    loss = "mean_squared_error",
    optimizer = "adam"
  )

  # fit model
  Sampling <- sample(nrow(x_train))
  model %>% keras::fit(
    x = x_train[Sampling[1:nbSubsamples.NLPCA], ],
    y = x_train[Sampling[1:nbSubsamples.NLPCA], ],
    epochs = 100,
    verbose = 0
  )
  # # evaluate the performance of the model
  # mse.ae2 <- keras::evaluate(model, x_train, x_train)
  # mse.ae2
  # extract the bottleneck layer
  intermediate_layer_model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output)
  intermediate_output <- predict(intermediate_layer_model, x_train)
  my_list <- list("Model" = model, "dataPCA" = intermediate_output, "MinVal" = Subset$MinVal, "MaxVal" = Subset$MaxVal)
  return(my_list)
}

# applies NLPCA to data
#
# @param Data matrix to apply NLPCA on
# @param model.AutoEncod model defined previously
#
# @return Image.Chunk after NLPCA
#' @importFrom future.apply future_lapply
predict.from.NLPCA <- function(Data, model.AutoEncod) {
  Image.Chunk <- future_lapply(Data, FUN = predict.from.NLPCA, model = model.AutoEncod, future.scheduling = 1.0)
  return(Image.Chunk)
}

# Function to perform PCA on a matrix
#
# @param X matrix to apply PCA on
# @param type PCA (no rescale) or SPCA (rescale)
#
# @return list of PCA parameters (PCs from X, mean, eigenvectors and values)
pca <- function(X, type) {
  p <- ncol(X)
  if (type == "SPCA") {
    modPCA <- prcomp(X, scale = TRUE)
    sig <- modPCA$scale
  } else if (type == "PCA") {
    modPCA <- prcomp(X, scale = FALSE)
    sig <- rep(1, p)
  }
  my_list <- list("dataPCA" = modPCA$x, "mu" = modPCA$center, "scale" = sig, "eiV" = modPCA$rotation, "eig" = modPCA$sdev)
  return(my_list)
}

# defines the number of pixels per iteration based on a trade-off between image size and sample size per iteration
#
# @param ImNames Path and name of the images to be processed
# @param NbIter number of iterations peformed to average diversity indices
#
# @return Pix.Per.Iter number of pixels per iteration
Define.Pixels.Per.Iter <- function(ImNames, NbIter = NbIter) {
  ImPath <- ImNames$Input.Image
  ImPathShade <- ImNames$Mask.list
  # define dimensions of the image
  ImPathHDR <- Get.HDR.Name(ImPath)
  HDR <- read.ENVI.header(ImPathHDR)
  Image.Format <- ENVI.Type2Bytes(HDR)
  ipix <- as.double(HDR$lines)
  jpix <- as.double(HDR$samples)
  Nb.Pixels <- ipix * jpix
  lenTot <- Nb.Pixels * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image.Format$Bytes) / (1024^3)
  # if shade mask, update number of pixels
  if ((!ImPathShade == FALSE) & (!ImPathShade == "")) {
    # read shade mask
    fid <- file(
      description = ImPathShade, open = "rb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    ShadeMask <- readBin(fid, integer(), n = Nb.Pixels, size = 1)
    close(fid)
    ShadeMask <- aperm(array(ShadeMask, dim = c(jpix, ipix)))
    # number of sunlit pixels
    Nb.Pixels.Sunlit <- length(which(ShadeMask == 1))
  } else {
    Nb.Pixels.Sunlit <- Nb.Pixels
  }

  # adjust the number of pixels per iteration
  # trade-off between number of pixels and total pixel size
  # maximum number of pixels to be used
  Max.Pixel.Per.Iter <- 250000
  Max.Pixel.Per.Iter.Size <- NbIter * (Max.Pixel.Per.Iter * as.double(HDR$bands) * Image.Format$Bytes) / (1024^3)
  # maximum amount of data (Gb) to be used
  Max.Size.Per.Iter <- 0.3
  # amount of data available after first mask
  ImDataGb <- (Nb.Pixels.Sunlit * as.double(HDR$bands) * Image.Format$Bytes) / (1024^3)

  # if Max.Pixel.Per.Iter correspond to reasonable size (<Max.Size.Per.Iter)
  if (Max.Pixel.Per.Iter.Size < Max.Size.Per.Iter) {
    # if enough data
    if (ImDataGb >= Max.Pixel.Per.Iter.Size) {
      Pix.Per.Iter <- Max.Pixel.Per.Iter
    } else if (ImDataGb < Max.Pixel.Per.Iter.Size) {
      # define number of pixels corresponding to ImDataGb
      Pix.Per.Iter <- floor(Nb.Pixels.Sunlit / NbIter)
    }
    # if size too important, adjust number of pixels to match Max.Size.Per.Iter
  } else if (Max.Pixel.Per.Iter.Size >= Max.Size.Per.Iter) {
    Pix.Per.Iter <- floor((Max.Size.Per.Iter * (1024^3) / (as.double(HDR$bands) * Image.Format$Bytes)) / NbIter)
  }
  return(Pix.Per.Iter)
}

#' Check if principal components are properly selected as expected by the method
#'
#' @param Output.Dir output directory
#' @param Input.Image.File path for image to be processed
#' @param PCA.Files path of PCA files
#' @param TypePCA Type of PCA (PCA, SPCA, NLPCA...)
#'
#' @return nothing
#' @export
Select.Components <- function(Input.Image.File, Output.Dir, PCA.Files, TypePCA = "SPCA", File.Open = FALSE) {
  message("")
  message("*********************************************************")
  message("Please check following PCA file:")
  print(PCA.Files)
  message("*********************************************************")
  Image.Name <- strsplit(basename(Input.Image.File), "\\.")[[1]][1]
  Output.Dir.Full <- paste(Output.Dir, "/", Image.Name, "/", TypePCA, "/", sep = "")
  Sel.PC <- paste(Output.Dir.Full, "/PCA/Selected_Components.txt", sep = "")
  message("list the principal components that will be used to estimate biodiversity in the file")
  message("")
  print(Sel.PC)
  message("")
  message("Then press ENTER")

  if (!file.exists(Sel.PC)) {
    file.create(Sel.PC)
  }
  if (File.Open == TRUE) {
    file.edit(Sel.PC, title=basename(Sel.PC))
  }
  message("*********************************************************")
  message("")
  readline(prompt = "")

  # if (!file.exists(Sel.PC)){
  #   print(Sel.PC)
  #   print('Does not exist')
  #   print('Please create tis file base on visual interpretation of')
  #   print(PCA.Files)
  #   print('Then press ENTER')
  #   readline(prompt = "")
  # }
  return()
}

# this function performs rescaling and
# either defines min and max from each feature in a data set,
# or applies the transformation based on a previously defined min and max
#
# @param x data matrix
# @param mode 'define' or 'apply'
# @param MinX  if 'apply'
# @param MaxX  if 'apply'
#
# @return rescaled data, min and max values
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
