# ==============================================================================
# biodivMapR
# Lib_ImageProcess.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
# ==============================================================================
# This Library contains functions to manipulate & process raster images
# Mainly applicable to ENVI HDR data wth BIL interleave
# ==============================================================================

# rebuild full image from list of subsets
#
# @param Image.Subsets subsets of an image. list expected
# @param ipix nb of lines for the full image
# @param jpix nb of columns for the full image
# @param nbBands nb of bands for the full image & subsets
#
# @return image full size in 3 dimensions
Build.Image.From.List <- function(Image.Subsets, ipix, jpix, nbBands) {
  Image <- array(NA, c(ipix, jpix, nbBands))
  Line.Begin <- 0
  Line.End <- 0
  for (i in 1:length(Image.Subsets)) {
    Line.Begin <- Line.End + 1
    Line.End <- Line.End + dim(Image.Subsets[[i]])[1]
    Image[Line.Begin:Line.End, , ] <- Image.Subsets[[i]]
  }
  rm(Image.Subsets)
  gc()
  return(Image)
}

# center and reduce data matrix based on known mean and SD
#
# @param X data matrix (each column is centered/reduced)
# @param m mean of each variable in the data matrix
# @param sig SD of each variable in the data matrix
#
# @return Centered matrix
Center.Reduce <- function(X, m, sig) {
  for (i in 1:ncol(X)) {
    X[, i] <- (X[, i] - m[i]) / sig[i]
  }
  return(X)
}

# change resolution in a HDR file
#
# @param HDR information read from a header file
# @param window_size multiplying factor for initial resolution
#
# @return updated HDR information
Change.Resolution.HDR <- function(HDR, window_size) {
  MapInfo <- strsplit(HDR$`map info`, split = ",")
  MapInfo[[1]][6] <- as.numeric(MapInfo[[1]][6]) * window_size
  MapInfo[[1]][7] <- as.numeric(MapInfo[[1]][7]) * window_size
  HDR$`map info` <- paste(MapInfo[[1]], collapse = ",")
  return(HDR)
}

# clean data matrix from inf values and identifies if some values
# (columns) are constant. variables showing constant value
# need to be eliminated before PCA
#
# @param DataMatrix each variable is a column
# @param Spectral summary of spectral information: which spectral bands selected from initial data
#
# @return updated DataMatrix and Spectral
Check.Data <- function(DataMatrix, Spectral) {
  # samples with inf value are eliminated
  for (i in 1:ncol(DataMatrix)) {
    elim <- which(DataMatrix[, i] == Inf)
    if (length(elim) > 0) {
      DataMatrix <- DataMatrix[-elim, ]
    }
  }
  # bands showing null std are eliminated
  stdsub <- apply(DataMatrix, 2, sd)
  BandsNoVar <- which(stdsub == 0)
  # BandsNoVar  = which(stdsub<=0.002)
  if (!length(Spectral$Bands2Keep[BandsNoVar]) == 0) {
    DataMatrix <- DataMatrix[, -BandsNoVar]
  }
  # !! the wl which is discarded correspond to original spectral bands,
  # whereas BandsNoVar corresponds to spectral band after contiuum removal
  Spectral$BandsNoVar <- BandsNoVar
  my_list <- list("DataMatrix" = DataMatrix, "Spectral" = Spectral)
  return(my_list)
}

# create a binary mask file based on a matrix and the corresponding header
#
# @param Mask matrix of a binary mask
# @param HDR header corresponding to an image to be masked
# @param Path.Mask path for the mask
#
# @return
Create.Shademask <- function(Mask, HDR, Path.Mask) {
  ipix <- HDR$lines
  jpix <- HDR$samples
  Mask <- array(Mask, c(ipix, jpix, 1))
  Mask <- aperm(Mask, c(2, 3, 1))
  fidOUT <- file(
    description = Path.Mask, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  writeBin(c(as.integer(Mask)), fidOUT, size = 1, endian = .Platform$endian, useBytes = FALSE)
  close(fidOUT)
  # write updated shademask
  HDR.Update <- HDR
  HDR.Update$bands <- 1
  HDR.Update$`data type` <- 1
  HDR.Update$`band names` <- {
    "Mask"
  }
  HDR.Update$wavelength <- NULL
  HDR.Update$fwhm <- NULL
  HDR.Update$resolution <- NULL
  HDR.Update$bandwidth <- NULL
  HDR.Update$purpose <- NULL
  HDRFpath <- paste(Path.Mask, ".hdr", sep = "")
  write.ENVI.header(HDR.Update, HDRFpath)
  return()
}

# define output directory and create it if necessary
#
# @param Output.Dir output directory
# @param ImPath image path
# @param TypePCA Type of PCA (PCA, SPCA, NLPCA...)
#
# @return path of the output directory
Define.Output.Dir <- function(Output.Dir, ImPath, TypePCA) {
  Image.Name <- strsplit(basename(ImPath), "\\.")[[1]][1]
  Output.Dir <- paste(Output.Dir, "/", Image.Name, "/", TypePCA, "/", sep = "")
  dir.create(Output.Dir, showWarnings = FALSE, recursive = TRUE)
  return(Output.Dir)
}

# define output directory and subdirectory and create it if necessary
#
# @param Output.Dir output directory
# @param ImPath image path
# @param TypePCA Type of PCA (PCA, SPCA, NLPCA...)
# @param Sub subdirectory
#
# @return path of the output directory
Define.Output.SubDir <- function(Output.Dir, ImPath, TypePCA, Sub) {
  Image.Name <- strsplit(basename(ImPath), "\\.")[[1]][1]
  Output.Dir <- paste(Output.Dir, "/", Image.Name, "/", TypePCA, "/", Sub, "/", sep = "")
  dir.create(Output.Dir, showWarnings = FALSE, recursive = TRUE)
  return(Output.Dir)
}

# get information corresponding to data type defined in ENVI
#
# @param Header Path of the hdr file
#
# @return description of data format corresponding to ENVI type
ENVI.Type2Bytes <- function(Header) {

  # http://www.harrisgeospatial.com/docs/ENVIHeaderFiles.html
  DataTypeImage <- Header$`data type`
  if (DataTypeImage == 1) {
    #   1 = Byte: 8-bit unsigned integer
    nbBytes <- 1
    Type <- "INT"
    is.Signed <- FALSE
  } else if (DataTypeImage == 2) {
    #   2 = Integer: 16-bit signed integer
    nbBytes <- 2
    Type <- "INT"
    is.Signed <- TRUE
  } else if (DataTypeImage == 3) {
    #   3 = Long: 32-bit signed integer
    nbBytes <- 4
    Type <- "INT"
    is.Signed <- TRUE
  } else if (DataTypeImage == 4) {
    #   4 = Floating-point: 32-bit single-precision
    nbBytes <- 4
    Type <- "FLOAT"
    is.Signed <- TRUE
  } else if (DataTypeImage == 5) {
    #   5 = Double-precision: 64-bit double-precision floating-point
    nbBytes <- 8
    Type <- "FLOAT"
    is.Signed <- TRUE
  } else if (DataTypeImage == 12) {
    #   12 = Unsigned integer: 16-bit
    nbBytes <- 2
    Type <- "INT"
    is.Signed <- FALSE
  } else if (DataTypeImage == 13) {
    #   13 = Unsigned long integer: 32-bit
    nbBytes <- 4
    Type <- "INT"
    is.Signed <- FALSE
  } else if (DataTypeImage == 14) {
    #   14 = 64-bit long integer (signed)
    nbBytes <- 8
    Type <- "INT"
    is.Signed <- TRUE
  } else if (DataTypeImage == 15) {
    #   15 = 64-bit unsigned long integer (unsigned)
    nbBytes <- 8
    Type <- "INT"
    is.Signed <- FALSE
  }
  if (Header$`byte order` == 0) {
    ByteOrder <- "little"
  } else if (Header$`byte order` == 1) {
    ByteOrder <- "big"
  }
  my_list <- list("Bytes" = nbBytes, "Type" = Type, "Signed" = is.Signed, "ByteOrder" = ByteOrder)
  return(my_list)
}

# define Water Vapor bands based on spectral smapling of original image
#
# @param ImPath path of the image
# @param Excluded.WL spectral domains corresponding to water vapor absorption
#
# @return bands corresponding to atmospheric water absorption domain
Exclude.Spectral.Domains <- function(ImPath, Excluded.WL = FALSE) {
  # definition of water vapor absorption
  if (length(Excluded.WL) == 1) {
    if (Excluded.WL == FALSE) {
      Excluded.WL <- c(0, 400)
      Excluded.WL <- rbind(Excluded.WL, c(895, 1005))
      Excluded.WL <- rbind(Excluded.WL, c(1180, 1480))
      Excluded.WL <- rbind(Excluded.WL, c(1780, 2040))
    }
  }
  # get image header data
  ImPathHDR <- Get.HDR.Name(ImPath)
  HDR <- read.ENVI.header(ImPathHDR)
  nbchannels0 <- HDR$bands
  idOkBand <- seq(1, nbchannels0)
  if ("wavelength" %in% names(HDR)) {
    wl <- HDR$wavelength
    WaterVapor <- c()
    for (w in 1:nrow(Excluded.WL)) {
      WaterVapor <- c(WaterVapor, which(wl > Excluded.WL[w, 1] & wl < Excluded.WL[w, 2]))
    }
    if (!length(WaterVapor) == 0) {
      wl <- wl[-WaterVapor]
      idOkBand <- idOkBand[-WaterVapor]
    }
  } else {
    wl <- integer()
    WaterVapor <- integer()
    idOkBand <- seq(1, nbchannels0)
  }
  my_list <- list("Wavelength" = wl, "WaterVapor" = WaterVapor, "Bands2Keep" = idOkBand)
  return(my_list)
}

# extracts sample points from binary image file
#
# @param coordPix.List coordinates of pixels to sample
# @param ImPath path for image
# @param HDR hdr path
#
# @return samples from image subset corresponding to coordPix.List
Extract.Pixels <- function(coordPix.List, ImPath, HDR) {
  coordPix.List <- matrix(coordPix.List, ncol = 2)
  Sample.Sel <- matrix(0, nrow = nrow(coordPix.List), ncol = HDR$bands)
  # each line of the image is read and the datapoints included in this line are extracted
  # seek file until the first line
  if (dim(coordPix.List)[1] > 1) {
    minRow <- min(coordPix.List[, 1])
    maxRow <- max(coordPix.List[, 1])
  } else if (dim(coordPix.List)[1] == 1) {
    minRow <- min(coordPix.List[1])
    maxRow <- max(coordPix.List[1])
  }
  nbRows <- maxRow - minRow + 1
  Pix.Per.Line <- as.double(HDR$samples) * as.double(HDR$bands)
  Pix.Per.Read <- as.double(nbRows) * as.double(HDR$samples) * as.double(HDR$bands)
  Image.Format <- ENVI.Type2Bytes(HDR)
  # open file and read section of interest
  fid <- file(
    description = ImPath, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!minRow == 1) {
    # nb of bits to skip
    nbSkip <- as.double(minRow - 1) * as.double(HDR$samples) * as.double(HDR$bands) * Image.Format$Bytes
    seek(fid, where = nbSkip, origin = "start", rw = "read")
  }
  if (Image.Format$Type == "INT") {
    ImgSubset <- readBin(fid, integer(), n = as.double(nbRows) * as.double(HDR$samples) * as.double(HDR$bands), size = Image.Format$Bytes, signed = Image.Format$Signed, endian = Image.Format$ByteOrder)
  } else if (Image.Format$Type == "FLOAT") {
    ImgSubset <- readBin(fid, numeric(), n = as.double(nbRows) * as.double(HDR$samples) * as.double(HDR$bands), size = Image.Format$Bytes, signed = Image.Format$Signed, endian = Image.Format$ByteOrder)
  }
  close(fid)
  # reshape ImgSubset in a 2D matrix in order to get selected pixels based on index
  ImgSubset <- array(aperm(array(ImgSubset, dim = c(HDR$samples, HDR$bands, nbRows)), c(3, 1, 2)), c(HDR$samples * nbRows, HDR$bands))
  # coordinates of the samples correspond to the whole image: need to convert to fit sample
  coordPix.List[, 1] <- coordPix.List[, 1] - minRow + 1
  # Get index of each pixel
  IndPix <- (coordPix.List[, 2] - 1) * nbRows + coordPix.List[, 1]
  # Get samples from subset
  Sample.Sel <- ImgSubset[IndPix, ]
  rm(ImgSubset)
  gc()
  return(Sample.Sel)
}

# extracts pixels from image based on their coordinates
#
# @param ImPath path for image
# @param MaxRAM
# @param ProgressBar
# @param Already.Split
# @param coordPix pixel coordinates
#
# @return samples from image corresponding to coordPix
Extract.Samples.From.Image <- function(ImPath, coordPix, MaxRAM = FALSE, ProgressBar = FALSE, Already.Split = FALSE) {
  # get image header data
  ImPathHDR <- Get.HDR.Name(ImPath)
  HDR <- read.ENVI.header(ImPathHDR)

  # compute the ranking of initial pixel list compared to index ranking
  if (typeof(coordPix) == "double" & dim(coordPix)[2] == 2) {
    if (dim(coordPix)[1] >= 2) {
      coordPix.tmp <- list()
      coordPix.tmp$Row <- coordPix[, 1]
      coordPix.tmp$Column <- coordPix[, 2]
    } else if (dim(coordPix)[1] == 1) {
      coordPix.tmp <- list()
      coordPix.tmp$Row <- coordPix[1]
      coordPix.tmp$Column <- coordPix[2]
    }
  } else if (typeof(coordPix) == "list" & length(grep("Row", names(coordPix))) > 0 & length(grep("Column", names(coordPix))) > 0) {
    coordPix.tmp <- coordPix
  }
  # initial index value of the pixels requested in the image, following original ranking
  Index.Init <- sub2ind(HDR, coordPix.tmp)
  # rank of the initial list of pixels
  Rank.Index.Init <- sort(Index.Init, index = TRUE)
  # convert
  if (typeof(coordPix) == "list" & length(grep("Row", names(coordPix))) > 0 & length(grep("Column", names(coordPix))) > 0) {
    coordPix <- cbind(coordPix$Row, coordPix$Column)
  }
  # divide data access based on the size of the image (to avoid reading 30 Gb file at once...)
  Image.Format <- ENVI.Type2Bytes(HDR)
  lenTot <- as.double(HDR$lines) * as.double(HDR$samples) * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image.Format$Bytes) / (1024^3)
  # maximum image size read at once. If image larger, then reads in multiple pieces
  if (MaxRAM == FALSE) {
    LimitSizeGb <- 0.25
  } else {
    LimitSizeGb <- MaxRAM
  }
  if (ImSizeGb < LimitSizeGb | Already.Split == TRUE) {
    Lines.Per.Read <- HDR$lines
    nbPieces <- 1
  } else {
    # nb of lines corresponding to LimitSizeGb
    OneLine <- as.double(HDR$samples) * as.double(HDR$bands) * Image.Format$Bytes
    Lines.Per.Read <- floor(LimitSizeGb * (1024^3) / OneLine)
    # number of pieces to split the image into
    nbPieces <- ceiling(HDR$lines / Lines.Per.Read)
  }

  # here split based on even number of pxiels to sample
  # should be based on image size in order to avoid memory
  # problems if target pixels  unevenly distributed in image
  coordPix.List <- Split.Pixel.Samples(coordPix, Lines.Per.Read)
  Sample.Sel <- list()

  # print('Extracting samples from image')
  # print(ImPath)
  if (ProgressBar == TRUE) {
    pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  }
  for (i in 1:length(coordPix.List)) {
    Sample.Sel[[i]] <- Extract.Pixels(coordPix.List[[i]], ImPath, HDR)
    if (ProgressBar == TRUE) {
      setTxtProgressBar(pb, 100 * i / length(coordPix.List))
    }
  }
  if (ProgressBar == TRUE) {
    close(pb)
  }
  Sample.Sel <- do.call("rbind", Sample.Sel)
  coordPix.List <- do.call("rbind", coordPix.List)

  # re-sort samples
  Sample.Sort <- list()
  Sample.Sort$Row <- coordPix.List[, 1]
  Sample.Sort$Column <- coordPix.List[, 2]
  Coord.Pixels <- sub2ind(HDR, Sample.Sort)
  # rank of the pixels extracted from the image
  Rank.Pixels <- sort(Coord.Pixels, index = TRUE)
  # pixel re-ordering needs to be performed in order to get back to Rank.Index.Init
  # first re-order pixels in order to follow increasing index value
  # Sample.Sort$index   = Coord.Pixels[Rank.Pixels$ix]
  # Sample.Sort$Row     = Sample.Sort$Row[Rank.Pixels$ix]
  # Sample.Sort$Column  = Sample.Sort$Column[Rank.Pixels$ix]
  # if bug, check this line
  # Sample.Sel          = Sample.Sel[Rank.Pixels$ix]
  Sample.Sel <- Sample.Sel[Rank.Pixels$ix, ]

  # then apply initial ranking as defined by Rank.Index.Init
  if (dim(Sample.Sel)[1] > 1) {
    Sample.Sel[Rank.Index.Init$ix, ] <- Sample.Sel
  } else if (dim(Sample.Sel)[1] == 1) {
    Sample.Sel <- Sample.Sel
  }
  return(Sample.Sel)
}




# Perform random sampling on an image including a shade mask
#
# @param ImPath path for image
# @param HDR path for hdr file
# @param ImPathShade path for shade mask
# @param NbIter number of iterations
# @param Pix.Per.Iter number of pixels per iteration
#
# @return samples from image and updated number of pixels to sampel if necessary
Get.Random.Subset.From.Image <- function(ImPath, HDR, ImPathShade, NbIter, Pix.Per.Iter) {
  nbPix2Sample <- NbIter * Pix.Per.Iter
  # get total number of pixels
  nbpix <- as.double(HDR$lines) * as.double(HDR$samples)
  # 1- Exclude masked pixels from random subset
  # Read Mask
  if ((!ImPathShade == "") & (!ImPathShade == FALSE)) {
    fid <- file(
      description = ImPathShade, open = "rb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    ShadeMask <- readBin(fid, integer(), n = nbpix, size = 1)
    close(fid)
    ShadeMask <- aperm(array(ShadeMask, dim = c(HDR$samples, HDR$lines)))
  } else {
    ShadeMask <- array(ones(HDR$lines * HDR$samples, 1), dim = c(HDR$lines, HDR$samples))
  }
  # get a list of samples among unmasked pixels
  IndexInit <- (matrix(1:nbpix, ncol = HDR$samples))
  ValidPixels <- sample(which(ShadeMask == 1 | ShadeMask == 255))
  NbValidPixels <- length(ValidPixels)
  # Check if number of valid pixels is compatible with number of pixels to be extracted
  # if number of pixels to sample superior to number of valid pixels, then adjust iterations
  if (NbValidPixels < nbPix2Sample) {
    nbPix2Sample <- NbValidPixels
    NbIter <- ceiling(NbValidPixels / Pix.Per.Iter)
    Pix.Per.Iter <- floor(NbValidPixels / NbIter)
    nbPix2Sample <- NbIter * Pix.Per.Iter
  }
  # Select a subset of nbPix2Sample among pixselected
  pixselected <- ValidPixels[1:nbPix2Sample]

  # define location of pixselected in binary file (avoid multiple reads) and optimize access to disk
  # the file is a BIL file, which means that for each line in the image,
  # we first need to define if datapoints are on the line, then read one point after the other
  coordi <- ((pixselected - 1) %% HDR$lines) + 1
  coordj <- floor((pixselected - 1) / HDR$lines) + 1
  # sort based on line and col
  indxPix <- order(coordi, coordj, na.last = FALSE)
  coordPix <- cbind(coordi[indxPix], coordj[indxPix])

  # 2- Extract samples from image
  Sample.Sel <- Extract.Samples.From.Image(ImPath, coordPix)
  # randomize!
  Sample.Sel <- Sample.Sel[sample(dim(Sample.Sel)[1]), ]
  my_list <- list("DataSubset" = Sample.Sel, "nbPix2Sample" = nbPix2Sample)
  return(my_list)
}

# does the system work with little endians or big endians?
#
# @return ByteOrder
#' @import tools
Get.Byte.Order <- function() {
  if (.Platform$endian == "little") {
    ByteOrder <- 0
  } else if (.Platform$endian == "big") {
    ByteOrder <- 1
  }
  return(ByteOrder)
}

# get hdr name from image file name, assuming it is BIL format
#
# @param ImPath path of the image
#
# @return corresponding hdr
#' @import tools
Get.HDR.Name <- function(ImPath) {
  if (file_ext(ImPath) == "") {
    ImPathHDR <- paste(ImPath, ".hdr", sep = "")
  } else if (file_ext(ImPath) == "bil") {
    ImPathHDR <- gsub(".bil", ".hdr", ImPath)
  } else if (file_ext(ImPath) == "zip") {
    ImPathHDR <- gsub(".zip", ".hdr", ImPath)
  } else {
    ImPathHDR <- paste(file_path_sans_ext(ImPath), ".hdr", sep = "")
  }

  if (!file.exists(ImPathHDR)) {
    message("WARNING : COULD NOT FIND HDR FILE")
    print(ImPathHDR)
    message("Process may stop")
  }
  return(ImPathHDR)
}

# gets rank of spectral bands in an image
#
# @param Spectral.Bands wavelength (nm) of the spectral bands to be found
# @param wavelength wavelength (nm) of all wavelengths in the image
#
# @return rank of all spectral bands of interest in the image and corresponding wavelength
Get.Image.Bands <- function(Spectral.Bands, wavelength) {
  ImBand <- c()
  Distance2WL <- c()
  for (band in Spectral.Bands) {
    Closest.Band <- order(abs(wavelength - band))[1]
    ImBand <- c(ImBand, Closest.Band)
    Distance2WL <- c(Distance2WL, abs(wavelength[Closest.Band] - band))
  }
  my_list <- list("ImBand" = ImBand, "Distance2WL" = Distance2WL)
  return(my_list)
}

# convert image coordinates from index to X-Y
#
# @param Raster image raster object
# @param Image.Index coordinates corresponding to the raster
ind2sub <- function(Raster, Image.Index) {
  c <- ((Image.Index - 1) %% Raster@ncols) + 1
  r <- floor((Image.Index - 1) / Raster@ncols) + 1
  my_list <- list("Column" = c, "Row" = r)
  return(my_list)
}

# convert image coordinates from index to X-Y
# image coordinates are given as index = (ID.col-1) * total.lines + ID.row
#
# @param Raster image raster object
# @param Image.Index coordinates corresponding to the raster
ind2sub2 <- function(Raster, Image.Index) {
  r <- ((Image.Index - 1) %% Raster@nrows) + 1
  c <- floor((Image.Index - 1) / Raster@nrows) + 1
  my_list <- list("Column" = c, "Row" = r)
  return(my_list)
}

# applies mean filter to an image
#
# @param ImageInit image defined as 2d matrix
# @param nbi number of lines in image
# @param nbj number of columns in image
# @param SizeFilt size of the window to compute mean filter
#
# @return rank of all spectral bands of interest in the image and corresponding wavelength
#' @importFrom matlab padarray
Mean.Filter <- function(ImageInit, nbi, nbj, SizeFilt) {
  E <- padarray(ImageInit, c(SizeFilt, SizeFilt), "symmetric", "both")
  ImageSmooth <- matrix(0, nrow = nbi, ncol = nbj)
  Mat2D <- MatSun <- matrix(0, nrow = ((2 * SizeFilt) + 1)^2, ncol = nbj)
  spl <- split(1:nbj, 1:((2 * SizeFilt) + 1))
  mid <- ceiling((((2 * SizeFilt) + 1)^2) / 2)
  for (i in (SizeFilt + 1):(nbi + SizeFilt)) {
    for (j in 1:((2 * SizeFilt) + 1)) {
      # create a 2D matrix
      Mat2D[, spl[[j]]] <- matrix(E[(i - SizeFilt):(i + SizeFilt), (spl[[j]][1]):(tail(spl[[j]], n = 1) + 2 * SizeFilt)], nrow = ((2 * SizeFilt) + 1)^2)
    }
    ImageSmooth[(i - SizeFilt), ] <- colMeans(Mat2D, na.rm = TRUE)
  }
  return(ImageSmooth)
}

# reads a subset from a binary image
#
# @param Byte.Start location of byte where to start reading in the image
# @param nbLines number of lines to read
# @param lenBin number of elements to read
# @param ImPath path for the image
# @param ImBand bands of interest
# @param jpix number of columns in the image
# @param nbChannels total number of channels in the image
# @param Image.Format type of data (INT/FLOAT)
#
# @return data corresponding to the subset in original 3D format
Read.Bin.Subset <- function(Byte.Start, nbLines, lenBin, ImPath, ImBand, jpix, nbChannels, Image.Format) {
  # number of bands to be kept
  nbSubset <- length(ImBand)
  # open file
  fid <- file(
    description = ImPath, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!Byte.Start == 1) {
    # skip the beginning of the file of not wanted
    seek(fid, where = as.double(Image.Format$Bytes * (Byte.Start - 1)), origin = "start", rw = "read")
  }
  if (Image.Format$Type == "INT") {
    linetmp <- readBin(fid, integer(), n = lenBin, size = Image.Format$Bytes, endian = Image.Format$ByteOrder)
  } else if (Image.Format$Type == "FLOAT") {
    linetmp <- readBin(fid, numeric(), n = lenBin, size = Image.Format$Bytes, endian = Image.Format$ByteOrder)
  }
  close(fid)
  # reshape data into original image subset
  linetmp <- aperm(array(linetmp, dim = c(jpix, nbChannels, nbLines)), c(3, 1, 2))
  linetmp <- array(linetmp[, , ImBand], dim = c(nbLines, jpix, length(ImBand)))
  if (nbLines == 1) {
    linetmp <- array(linetmp, c(nbLines, jpix, nbSubset))
  }
  return(linetmp)
}

# Reads ENVI hdr file
#
# @param headerFpath Path of the hdr file
#
# @return list of the content of the hdr file
read.ENVI.header <- function(headerFpath) {
  # header <- paste(header, collapse = "\n")
  if (!grepl(".hdr$", headerFpath)) {
    stop("File extension should be .hdr")
  }
  header <- readLines(headerFpath)
  ## check ENVI at beginning of file
  if (!grepl("ENVI", header[1])) {
    stop("Not an ENVI header (ENVI keyword missing)")
  } else {
    header <- header [-1]
  }
  ## remove curly braces and put multi-line key-value-pairs into one line
  header <- gsub("\\{([^}]*)\\}", "\\1", header)
  l <- grep("\\{", header)
  r <- grep("\\}", header)

  if (length(l) != length(r)) {
    stop("Error matching curly braces in header (differing numbers).")
  }

  if (any(r <= l)) {
    stop("Mismatch of curly braces in header.")
  }

  header[l] <- sub("\\{", "", header[l])
  header[r] <- sub("\\}", "", header[r])

  for (i in rev(seq_along(l))) {
    header <- c(
      header [seq_len(l [i] - 1)],
      paste(header [l [i]:r [i]], collapse = "\n"),
      header [-seq_len(r [i])]
    )
  }

  ## split key = value constructs into list with keys as names
  header <- sapply(header, split.line, "=", USE.NAMES = FALSE)
  names(header) <- tolower(names(header))

  ## process numeric values
  tmp <- names(header) %in% c(
    "samples", "lines", "bands", "header offset", "data type",
    "byte order", "default bands", "data ignore value",
    "wavelength", "fwhm", "data gain values"
  )
  header [tmp] <- lapply(header [tmp], function(x) {
    as.numeric(unlist(strsplit(x, ",")))
  })

  return(header)
}

# read specific image bands from image
#
# @param ImPath Path of the image to read
# @param Header Header for the image
# @param ImBand Bands to be read
#
# @return Image.Subset information corresponding to ImBand
Read.Image.Bands <- function(ImPath, Header, ImBand) {
  # first get image format
  Image.Format <- ENVI.Type2Bytes(Header)
  ipix <- Header$lines
  jpix <- Header$samples
  nbChannels <- Header$bands
  nbSubset <- length(ImBand)
  # then open and read image
  # depending on image size, need to read in one or multiple times
  lenTot <- as.double(ipix) * as.double(jpix) * as.double(nbChannels)
  lenSubset <- as.double(ipix) * as.double(jpix) * as.double(nbSubset)
  ImSizeGb <- (lenTot * Image.Format$Bytes) / (1024^3)
  # maximum image size read at once. If image larger, then reads in multiple pieces
  LimitSizeGb <- 0.25
  if (ImSizeGb < LimitSizeGb) {
    nbLinesPerCPU <- ceiling(ipix)
    nbPieces <- 1
  } else {
    # nb of lines corresponding to LimitSizeGb
    OneLine <- as.double(jpix) * as.double(nbChannels) * Image.Format$Bytes
    nbLinesPerCPU <- floor(LimitSizeGb * (1024^3) / OneLine)
    # number of pieces to split the image into
    nbPieces <- ceiling(ipix / nbLinesPerCPU)
  }
  # Define segments of image to be read
  SeqRead.Image <- Where.To.Read(Header, nbPieces)
  # Read segments (subsets) of image
  Image.Subsets <- list()
  for (i in 1:nbPieces) {
    # number of elements to be read
    Byte.Start <- SeqRead.Image$ReadByte.Start[i]
    nbLines <- SeqRead.Image$Lines.Per.Chunk[i]
    lenBin <- SeqRead.Image$ReadByte.End[i] - SeqRead.Image$ReadByte.Start[i] + 1
    Image.Subsets[[i]] <- Read.Bin.Subset(Byte.Start, nbLines, lenBin, ImPath, ImBand, jpix, nbChannels, Image.Format)
  }
  # reshape image with original size and selected bands
  Image.Subset <- Build.Image.From.List(Image.Subsets, ipix, jpix, nbSubset)
  rm(Image.Subsets)
  gc()
  return(Image.Subset)
}

# reads subset of an image
#
# @param ImPath path for the image
# @param HDR header information corresponding to the image
# @param Byte.Start location of byte where to start reading in the image
# @param lenBin number of elements to read
# @param nbLines number of lines to read
# @param Image.Format type of data (INT/FLOAT)
# @param ImgFormat should output be 2D or 3D (original image format)?
#
# @return data corresponding to the subset in original 3D format
Read.Image.Subset <- function(ImPath, HDR, Byte.Start, lenBin, nbLines, Image.Format, ImgFormat) {
  fidIm <- file(
    description = ImPath, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )

  # read relevant portion
  if (!Byte.Start == 1) {
    # skip the beginning of the file of not wanted
    seek(fidIm, where = Image.Format$Bytes * (Byte.Start - 1), origin = "start", rw = "read")
  }
  if (Image.Format$Type == "INT") {
    if (HDR$`data type` == 1) {
      Image.Chunk <- readBin(fidIm, integer(), n = lenBin, size = Image.Format$Bytes, signed = FALSE, endian = Image.Format$ByteOrder)
    } else {
      Image.Chunk <- readBin(fidIm, integer(), n = lenBin, size = Image.Format$Bytes, signed = TRUE, endian = Image.Format$ByteOrder)
    }
  } else if (Image.Format$Type == "FLOAT") {
    Image.Chunk <- readBin(fidIm, numeric(), n = lenBin, size = Image.Format$Bytes, endian = Image.Format$ByteOrder)
  }

  if (ImgFormat == "3D" | ImgFormat == "2D") {
    Image.Chunk <- aperm(array(Image.Chunk, dim = c(HDR$samples, HDR$bands, nbLines)), c(3, 1, 2))
  }
  if (ImgFormat == "2D") {
    Image.Chunk <- array(Image.Chunk, c(nbLines * HDR$samples, HDR$bands))
  }
  if (ImgFormat == "Shade") {
    Image.Chunk <- matrix(Image.Chunk, nbLines, HDR$samples, byrow = T)
  }
  close(fidIm)
  return(Image.Chunk)
}

# This image reads information in a zipfile
#
# @param Zip.Path path for the zipped binary image file in .BIL format and
# @param HDR Header corresponding
# @return Img.Data vecotr containing all image data
Read.Image.Zipfile <- function(Zip.Path) {

  # get HDR corresponding to zipfile
  # ptm <- proc.time()
  ImPathHDR <- Get.HDR.Name(Zip.Path)
  HDR <- read.ENVI.header(ImPathHDR)
  nb.Elements <- HDR$samples * HDR$lines * HDR$bands
  Image.Format <- ENVI.Type2Bytes(HDR$`data type`)
  fid <- unz(description = Zip.Path, filename = zip_list(Zip.Path)$filename, open = "rb", encoding = getOption("encoding"))
  if (Image.Format$Type == "INT") {
    Img.Data <- readBin(fid, integer(), n = nb.Elements, size = Image.Format$Bytes, endian = "little")
  } else if (Image.Format$Type == "FLOAT") {
    Img.Data <- readBin(fid, numeric(), n = nb.Elements, size = Image.Format$Bytes, endian = "little")
  }
  Img.Data <- aperm(array(Img.Data, dim = c(HDR$samples, HDR$bands, HDR$lines)), c(3, 1, 2))
  Img.Data <- c(array(Img.Data[, , HDR$bands], dim = c(HDR$lines, HDR$samples, length(HDR$bands))))
  close(fid)
  # proc.time() - ptm
  return(Img.Data)
}



#' ENVI functions
#'
#' based on https://github.com/cran/hyperSpec/blob/master/R/read.ENVI.R
#' added wavelength, fwhm, ... to header reading
#' Title
#'
#' @param x character.
#' @param separator character
#' @param trim.blank boolean.
#'
#' @return list.
#' @export
split.line <- function(x, separator, trim.blank = TRUE) {
  tmp <- regexpr(separator, x)
  key <- substr(x, 1, tmp - 1)
  value <- substr(x, tmp + 1, nchar(x))
  if (trim.blank) {
    blank.pattern <- "^[[:blank:]]*([^[:blank:]]+.*[^[:blank:]]+)[[:blank:]]*$"
    key <- sub(blank.pattern, "\\1", key)
    value <- sub(blank.pattern, "\\1", value)
  }
  value <- as.list(value)
  names(value) <- key
  return(value)
}

# splits a set of pixels to be sampled in an image based on number of lines, not number of samples
#
# @param coordPix coordinates of pixels to sample
# @param Lines.Per.Read max number of lines per read for memory concerns
#
# @return coordPix.List list of pixel coordinates
Split.Pixel.Samples <- function(coordPix, Lines.Per.Read) {
  # maximum Lines.Per.Read
  if (dim(coordPix)[1] > 1) {
    nb.Lines <- max(coordPix[, 1]) - min(coordPix[, 1]) + 1
  } else if (dim(coordPix)[1] == 1) {
    nb.Lines <- max(coordPix[1]) - min(coordPix[1]) + 1
  }
  nb.Pieces <- ceiling(nb.Lines / Lines.Per.Read)

  if (dim(coordPix)[1] > 1) {
    Min.Line <- ceiling(seq(min(coordPix[, 1]), max(coordPix[, 1]), length.out = nb.Pieces + 1))
  } else if (dim(coordPix)[1] == 1) {
    Min.Line <- ceiling(seq(min(coordPix[1]), max(coordPix[1]), length.out = nb.Pieces + 1))
  }

  Max.Line <- Min.Line - 1
  Min.Line <- Min.Line[-length(Min.Line)]
  Max.Line <- Max.Line[-1]
  Max.Line[length(Max.Line)] <- Max.Line[length(Max.Line)] + 1
  coordPix.List <- list()
  ii <- 0
  for (i in 1:length(Min.Line)) {
    if (dim(coordPix)[1] > 1) {
      selpix <- which(coordPix[, 1] >= Min.Line[i] & coordPix[, 1] <= Max.Line[i])
    } else if (dim(coordPix)[1] == 1) {
      selpix <- which(coordPix[1] >= Min.Line[i] & coordPix[1] <= Max.Line[i])
    }
    if (length(selpix) > 1) {
      ii <- ii + 1
      coordPix.List[[ii]] <- coordPix[selpix, ]
    } else if (length(selpix) == 1) {
      ii <- ii + 1
      coordPix.List[[ii]] <- coordPix
    }
  }
  return(coordPix.List)
}

# updates an existing shade mask
#
# @param ImPathShade original shade mask (may not exist)
# @param Header header correpondingproviding general info about data format
# @param Mask data to be used in the mask
# @param ImPathShade.Update path for teh updated shade mask
#
# @return ImPathShade.Update
Update.Shademask <- function(ImPathShade, Header, Mask, ImPathShade.Update) {
  ipix <- Header$lines
  jpix <- Header$samples
  nbpix <- as.double(ipix) * as.double(jpix)
  # if ImPathShade provided
  if ((!ImPathShade == "") & (!ImPathShade == FALSE)) {
    # read shade mask
    fid <- file(
      description = ImPathShade, open = "rb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    lenBin <- nbpix
    ShadeMask0 <- readBin(fid, integer(), n = lenBin, size = 1)
    close(fid)
    ShadeMask0 <- aperm(array(ShadeMask0, dim = c(jpix, ipix)))
    # multiply by Mask
    Mask <- ShadeMask0 * Mask
  }
  Mask <- array(Mask, c(ipix, jpix, 1))
  Mask <- aperm(Mask, c(2, 3, 1))
  fidOUT <- file(
    description = ImPathShade.Update, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  writeBin(c(as.integer(Mask)), fidOUT, size = 1, endian = .Platform$endian, useBytes = FALSE)
  close(fidOUT)
  # write updated shademask
  Header.Update <- Header
  Header.Update$bands <- 1
  Header.Update$`data type` <- 1
  Header.Update$`band names` <- {
    "Mask"
  }
  Header.Update$wavelength <- NULL
  Header.Update$fwhm <- NULL
  Header.Update$resolution <- NULL
  Header.Update$bandwidth <- NULL
  Header.Update$purpose <- NULL
  headerFpath <- paste(ImPathShade.Update, ".hdr", sep = "")
  write.ENVI.header(Header.Update, headerFpath)
  return(ImPathShade.Update)
}

# defines which byte should be read for each part of an image split in nbPieces
#
# @param HDR header info
# @param nbPieces number of pieces resulting from image split
#
# @return location of the bytes corresponding to beginning and end of each piece, and corresponding number of lines
Where.To.Read <- function(HDR, nbPieces) {
  Data.Per.Line <- as.double(HDR$samples) * as.double(HDR$bands)
  lenTot <- as.double(HDR$lines) * Data.Per.Line
  # starting line for each chunk
  Start.Per.Chunk <- ceiling(seq(1, (HDR$lines + 1), length.out = nbPieces + 1))
  # elements in input data
  lb <- 1 + ((Start.Per.Chunk - 1) * Data.Per.Line)
  ub <- lb - 1
  ReadByte.Start <- lb[1:nbPieces]
  ReadByte.End <- ub[2:(nbPieces + 1)]
  Lines.Per.Chunk <- diff(Start.Per.Chunk)
  my_list <- list("ReadByte.Start" = ReadByte.Start, "ReadByte.End" = ReadByte.End, "Lines.Per.Chunk" = Lines.Per.Chunk)
  return(my_list)
}

# defines which byte should be read for each part of an image split in nbPieces
#
# @param HDR header info
# @param nbPieces number of pieces resulting from image split
# @return location of the bytes corresponding to beginning and end of each piece, and corresponding number of lines
Where.To.Read.Kernel <- function(HDR, nbPieces, SE.Size) {
  Data.Per.Line <- as.double(HDR$samples) * as.double(HDR$bands)
  lenTot <- as.double(HDR$lines) * Data.Per.Line
  # starting line for each chunk
  Start.Per.Chunk <- ceiling(seq(1, (HDR$lines + 1), length.out = nbPieces + 1))
  Start.Per.Chunk <- Start.Per.Chunk - Start.Per.Chunk %% SE.Size + 1
  # elements in input data
  lb <- 1 + ((Start.Per.Chunk - 1) * Data.Per.Line)
  ub <- lb - 1
  ReadByte.Start <- lb[1:nbPieces]
  ReadByte.End <- ub[2:(nbPieces + 1)]
  Lines.Per.Chunk <- diff(Start.Per.Chunk)
  my_list <- list("ReadByte.Start" = ReadByte.Start, "ReadByte.End" = ReadByte.End, "Lines.Per.Chunk" = Lines.Per.Chunk)
  return(my_list)
}

# defines which byte should be written for each part of an image split in nbPieces
#
# @param HDR.SS header info for SpectralSpecies file
# @param HDR.SSD header info for SpectralSpecies_Distribution file
# @param nbPieces number of pieces resulting from image split
# @param SE.Size
#
# @return location of the bytes corresponding to beginning and end of each piece, and corresponding number of lines
Where.To.Write.Kernel <- function(HDR.SS, HDR.SSD, nbPieces, SE.Size) {
  Data.Per.Line.SS <- as.double(HDR.SS$samples) * as.double(HDR.SS$bands)
  Data.Per.Line.SSD <- as.double(HDR.SSD$samples) * as.double(HDR.SSD$bands)

  # starting line for each chunk of spectral species
  Start.Per.Chunk <- ceiling(seq(1, (HDR.SS$lines + 1), length.out = nbPieces + 1))
  Start.Per.Chunk <- Start.Per.Chunk - Start.Per.Chunk %% SE.Size
  Start.Per.Chunk.SSD <- (Start.Per.Chunk / SE.Size) + 1

  # elements in input data
  Image.Format <- ENVI.Type2Bytes(HDR.SSD)
  lb_SSD <- 1 + (((Start.Per.Chunk.SSD - 1) * Data.Per.Line.SSD) * Image.Format$Bytes)
  ub_SSD <- lb_SSD - 1
  ReadByte.Start.SSD <- lb_SSD[1:nbPieces]
  ReadByte.End.SSD <- ub_SSD[2:(nbPieces + 1)]
  Lines.Per.Chunk.SSD <- diff(Start.Per.Chunk.SSD)

  my_list <- list("ReadByte.Start" = ReadByte.Start.SSD, "ReadByte.End" = ReadByte.End.SSD, "Lines.Per.Chunk" = Lines.Per.Chunk.SSD)
  return(my_list)
}

# writes ENVI hdr file
#
# @param header content to be written
# @param headerFpath Path of the hdr file
#
# @return
#' @importFrom stringr str_count
write.ENVI.header <- function(header, headerFpath) {
  h <- lapply(header, function(x) {
    if (length(x) > 1 || (is.character(x) && str_count(x, "\\w+") > 1)) {
      x <- paste0("{", paste(x, collapse = ","), "}")
    }
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c("ENVI", paste(names(header), h, sep = " = ")), con = headerFpath)
}

# convert image coordinates from X-Y to index
#
# @param HDR.Raster
# @param Pixels coordinates corresponding to the raster
#
# @return Image.Index
sub2ind <- function(HDR.Raster, Pixels) {
  Image.Index <- (Pixels$Column - 1) * HDR.Raster$lines + Pixels$Row
  return(Image.Index)
}

# defines the number of pieces resulting from image split
#
# @param HDR information extracted from a header
# @param LimitSizeGb maximum size of individual pieces of an image (in Gb)
#
# @return nbPieces number of pieces
Split.Image <- function(HDR, LimitSizeGb = FALSE) {
  Image.Format <- ENVI.Type2Bytes(HDR)
  lenTot <- as.double(HDR$samples) * as.double(HDR$lines) * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image.Format$Bytes) / (1024^3)
  # maximum image size read at once. If image larger, then reads in multiple pieces
  if (LimitSizeGb == FALSE) {
    LimitSizeGb <- 0.25
  }
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
  return(nbPieces)
}

# revert resolution in a HDR file
#
# @param HDR information read from a header file
# @param window_size multiplying factor for initial resolution
#
# @return updated HDR information
Revert.Resolution.HDR <- function(HDR, window_size) {
  MapInfo <- strsplit(HDR$`map info`, split = ",")
  MapInfo[[1]][6] <- as.numeric(MapInfo[[1]][6]) / window_size
  MapInfo[[1]][7] <- as.numeric(MapInfo[[1]][7]) / window_size
  HDR$`map info` <- paste(MapInfo[[1]], collapse = ",")
  return(HDR)
}

# Zips an image file
#
# @param ImagePath path for the image
ZipFile <- function(ImagePath) {
  PathRoot <- getwd()
  ImageDir <- dirname(ImagePath)
  ImageName <- basename(ImagePath)
  setwd(ImageDir)
  zip::zip(zipfile = paste0(ImageName, ".zip"), files = ImageName)
  file.remove(ImageName)
  setwd(PathRoot)
  return()
}
