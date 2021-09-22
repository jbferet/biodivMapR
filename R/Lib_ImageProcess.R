# ==============================================================================
# biodivMapR
# Lib_ImageProcess.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library contains functions to manipulate & process raster images
# Mainly applicable to ENVI HDR data wth BIL interleave
# ==============================================================================

#' rebuilds full image from list of subsets
#'
#' @param Image_Subsets list. subsets of an image
#' @param ipix numeric. nb of lines for the full image
#' @param jpix numeric. nb of columns for the full image
#' @param nbBands numeric. nb of bands for the full image & subsets
#'
#' @return Image numeric. full size in 3 dimensions
#' @export

build_image_from_list <- function(Image_Subsets, ipix, jpix, nbBands) {
  Image <- array(NA, c(ipix, jpix, nbBands))
  Line_Begin <- 0
  Line_End <- 0
  for (i in 1:length(Image_Subsets)) {
    Line_Begin <- Line_End + 1
    Line_End <- Line_End + dim(Image_Subsets[[i]])[1]
    Image[Line_Begin:Line_End, , ] <- Image_Subsets[[i]]
  }
  rm(Image_Subsets)
  gc()
  return(Image)
}

#' center and reduce data matrix based on known mean and SD
#'
#' @param X numeric. data matrix (each column is centered/reduced)
#' @param m numeric. mean of each variable in the data matrix
#' @param sig numeric. SD of each variable in the data matrix
#'
#' @return X numeric. Centered matrix
#' @export

center_reduce <- function(X, m, sig) {
  for (i in 1:ncol(X)) {
    X[, i] <- (X[, i] - m[i]) / sig[i]
  }
  return(X)
}

#' change resolution in a HDR file
#'
#' @param HDR information read from a header file
#' @param window_size multiplying factor for initial resolution
#'
#' @return updated HDR information
#' @export

change_resolution_HDR <- function(HDR, window_size) {
  MapInfo <- strsplit(HDR$`map info`, split = ",")
  MapInfo[[1]][6] <- as.numeric(MapInfo[[1]][6]) * window_size
  MapInfo[[1]][7] <- as.numeric(MapInfo[[1]][7]) * window_size
  HDR$`map info` <- paste(MapInfo[[1]], collapse = ",")
  return(HDR)
}

#' create a hdr file for a raster
#'
#' @param ImPath character. path of the raster image
#' @param Sensor character. name of the sensor
#' @param SpectralBands numeric. vector of spectral bands for the raster
#' @param BandName character. name for each band in the raster
#' @param WLunits character. wavelengths unit. 'Micrometers' or 'Nanometers'
#'
#' @return corresponding hdr
#' @importFrom raster raster hdr
#' @export

create_hdr <- function(ImPath, Sensor = 'unknown', SpectralBands = NULL,
                       BandName = NULL, WLunits = NULL) {

  # hdr file corresponding to the file
  ImPathHDR <- get_HDR_name(ImPath,showWarnings = FALSE)
  if (file.exists(ImPathHDR)) {
    message("WARNING : HDR FILE ALREADY EXISTS")
    print(ImPathHDR)
    message("It will be overwritten")
  }
  # create HDR
  r <- raster::raster(ImPath)
  raster::hdr(x = r,format = 'ENVI')
  # add bands
  HDR_input <- read_ENVI_header(ImPathHDR)
  HDR_Temp_Path <- system.file("extdata", "HDR", paste0(Sensor, ".hdr"), package = "biodivMapR")
  if (file.exists(HDR_Temp_Path)) {
    # get raster template corresponding to the sensor
    HDR_Template <- read_ENVI_header(HDR_Temp_Path)
    # Define spectral bands (central wavelength and band name)
    if (is.null(SpectralBands)){
      HDR_input$wavelength <- HDR_Template$wavelength
    } else {
      HDR_input$wavelength <- SpectralBands
    }
    if (is.null(BandName)){
      HDR_input$`band names` <- HDR_Template$`band names`
    } else {
      HDR_input$`band names` <- BandName
    }
    # Define sensor type
    HDR_input$`sensor type` <- Sensor
    if (is.null(WLunits)) {
      HDR_input$`wavelength units` <- HDR_Template$`wavelength units`
    } else {
      HDR_input$`wavelength units` <- WLunits
    }
  } else if (!file.exists(HDR_Temp_Path)){
    message(paste('No template exist for the sensor',Sensor,sep = ' '))
    HDR_input$wavelength <- SpectralBands
    HDR_input$`band names` <- BandName
    HDR_input$`sensor type` <- Sensor
    HDR_input$`wavelength units` <- WLunits
    if (is.null(SpectralBands)){
      message('Please provide sensor information to update header')
      message('Central wavelength for each spectral band expected in SpectralBands')

    }
  }
  # define visual stretch in the VIS domain
  if (HDR_input$`data type`==2 | HDR_input$`data type`==4){
    HDR_input$`default stretch` <- "0 1000 linear"
  }
  # write corresponding hdr file
  write_ENVI_header(HDR_input, ImPathHDR)
  return(ImPathHDR)
}


#' remove constant bands
#'
#' @param DataMatrix numeric. each variable is a column
#' @param Spectral list. summary of spectral information: which spectral bands selected from initial data
#'
#' @return updated DataMatrix and Spectral
#' @importFrom stats sd
#' @export

rm_invariant_bands <- function(DataMatrix, Spectral) {
  # samples with inf value are eliminated
  for (i in 1:ncol(DataMatrix)) {
    elim <- which(DataMatrix[, i] == Inf)
    if (length(elim) > 0) {
      DataMatrix <- DataMatrix[-elim, ]
    }
  }
  # bands showing null std are removed
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

#' define output directory and create it if necessary
#'
#' @param Output_Dir character. output directory
#' @param ImPath character. image path
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...)
#'
#' @return path of the output directory
#' @export

define_output_directory <- function(Output_Dir, ImPath, TypePCA) {
  Image_Name <- strsplit(basename(ImPath), "\\.")[[1]][1]
  Output_Dir <- file.path(Output_Dir, Image_Name, TypePCA)
  dir.create(Output_Dir, showWarnings = FALSE, recursive = TRUE)
  return(Output_Dir)
}

#' define output directory and subdirectory and create it if necessary
#'
#' @param Output_Dir character. output directory
#' @param ImPath character. image path
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...)
#' @param Sub character. subdirectory
#'
#' @return path of the output directory
#' @export

define_output_subdir <- function(Output_Dir, ImPath, TypePCA, Sub) {
  if (ImPath==FALSE){
    message('Please provide input image file')
    stop()
  }
  Image_Name <- strsplit(basename(ImPath), "\\.")[[1]][1]
  Output_Dir <- file.path(Output_Dir, Image_Name, TypePCA, Sub)
  dir.create(Output_Dir, showWarnings = FALSE, recursive = TRUE)
  return(Output_Dir)
}

#' get information corresponding to data type defined in ENVI
#'
#' @param HDR list. header file
#'
#' @return description of data format corresponding to ENVI type
#' @export

ENVI_type2bytes <- function(HDR) {

  # http://www.harrisgeospatial.com/docs/ENVIHeaderFiles.html
  DataTypeImage <- HDR$`data type`
  if (DataTypeImage == 1) {
    #   1 = Byte: 8-bit unsigned integer
    nbBytes <- 1
    Type <- "INT"
    is_Signed <- FALSE
  } else if (DataTypeImage == 2) {
    #   2 = Integer: 16-bit signed integer
    nbBytes <- 2
    Type <- "INT"
    is_Signed <- TRUE
  } else if (DataTypeImage == 3) {
    #   3 = Long: 32-bit signed integer
    nbBytes <- 4
    Type <- "INT"
    is_Signed <- TRUE
  } else if (DataTypeImage == 4) {
    #   4 = Floating-point: 32-bit single-precision
    nbBytes <- 4
    Type <- "FLOAT"
    is_Signed <- TRUE
  } else if (DataTypeImage == 5) {
    #   5 = Double-precision: 64-bit double-precision floating-point
    nbBytes <- 8
    Type <- "FLOAT"
    is_Signed <- TRUE
  } else if (DataTypeImage == 12) {
    #   12 = Unsigned integer: 16-bit
    nbBytes <- 2
    Type <- "INT"
    is_Signed <- FALSE
  } else if (DataTypeImage == 13) {
    #   13 = Unsigned long integer: 32-bit
    nbBytes <- 4
    Type <- "INT"
    is_Signed <- FALSE
  } else if (DataTypeImage == 14) {
    #   14 = 64-bit long integer (signed)
    nbBytes <- 8
    Type <- "INT"
    is_Signed <- TRUE
  } else if (DataTypeImage == 15) {
    #   15 = 64-bit unsigned long integer (unsigned)
    nbBytes <- 8
    Type <- "INT"
    is_Signed <- FALSE
  }
  if (HDR$`byte order` == 0) {
    ByteOrder <- "little"
  } else if (HDR$`byte order` == 1) {
    ByteOrder <- "big"
  }
  my_list <- list("Bytes" = nbBytes, "Type" = Type, "Signed" = is_Signed, "ByteOrder" = ByteOrder)
  return(my_list)
}

# define Water Vapor bands based on spectral smapling of original image
#
# @param ImPath path of the image
# @param Excluded_WL spectral domains corresponding to water vapor absorption
#
# @return bands corresponding to atmospheric water absorption domain
exclude_spectral_domains <- function(ImPath, Excluded_WL = FALSE) {
  # definition of water vapor absorption
  if (is.null(Excluded_WL)){
    Excluded_WL <- c(0, 0)
  } else if (length(Excluded_WL) == 1) {
    if (Excluded_WL == FALSE) {
      Excluded_WL <- c(0, 400)
      Excluded_WL <- rbind(Excluded_WL, c(895, 1005))
      Excluded_WL <- rbind(Excluded_WL, c(1320, 1480))
      Excluded_WL <- rbind(Excluded_WL, c(1780, 2040))
    }
  }
  message("Exclude spectral domains corresponding to water vapor absorption")
  message("The following spectral domains (min and max spectral band in nanometers) will be discarded")
  print(Excluded_WL)
  message("Please define the input variable Excluded_WL when calling perform_PCA")
  message("if you want to modify these spectral domains")

  # in case a unique specrtal domain is provided as excluded domain
  if (length(Excluded_WL)==2){
    Excluded_WL = matrix(Excluded_WL,ncol = 2)
  }
  # get image header data
  ImPathHDR <- get_HDR_name(ImPath)
  HDR <- read_ENVI_header(ImPathHDR)
  if (!is.null(HDR$`wavelength units`)){
    if ((HDR$`wavelength units`=='Micrometers') | (HDR$`wavelength units`=='micrometers')){
      Excluded_WL <- 0.001*Excluded_WL
    } else if (max(HDR$wavelength)<100){
      Excluded_WL <- 0.001*Excluded_WL
    }
  }
  nbchannels0 <- HDR$bands
  idOkBand <- seq(1, nbchannels0)
  if ("wavelength" %in% names(HDR)) {
    wl <- HDR$wavelength
    WaterVapor <- c()
    for (w in 1:nrow(Excluded_WL)) {
      WaterVapor <- c(WaterVapor, which(wl > Excluded_WL[w, 1] & wl < Excluded_WL[w, 2]))
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
# @param coordPix_List coordinates of pixels to sample
# @param ImPath path for image
# @param HDR hdr path
#
# @return samples from image subset corresponding to coordPix_List
extract_pixels <- function(coordPix_List, ImPath, HDR) {
  coordPix_List <- matrix(coordPix_List, ncol = 2)
  Sample_Sel <- matrix(0, nrow = nrow(coordPix_List), ncol = HDR$bands)
  # each line of the image is read and the datapoints included in this line are extracted
  # seek file until the first line
  if (dim(coordPix_List)[1] > 1) {
    minRow <- min(coordPix_List[, 1])
    maxRow <- max(coordPix_List[, 1])
  } else if (dim(coordPix_List)[1] == 1) {
    minRow <- min(coordPix_List[1])
    maxRow <- max(coordPix_List[1])
  }
  nbRows <- maxRow - minRow + 1
  Pix_Per_Line <- as.double(HDR$samples) * as.double(HDR$bands)
  Pix_Per_Read <- as.double(nbRows) * as.double(HDR$samples) * as.double(HDR$bands)
  Image_Format <- ENVI_type2bytes(HDR)
  # open file and read section of interest
  fid <- file(
    description = ImPath, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!minRow == 1) {
    # nb of bits to skip
    nbSkip <- as.double(minRow - 1) * as.double(HDR$samples) * as.double(HDR$bands) * Image_Format$Bytes
    seek(fid, where = nbSkip, origin = "start", rw = "read")
  }
  if (Image_Format$Type == "INT") {
    ImgSubset <- readBin(fid, integer(), n = as.double(nbRows) * as.double(HDR$samples) * as.double(HDR$bands), size = Image_Format$Bytes, signed = Image_Format$Signed, endian = Image_Format$ByteOrder)
  } else if (Image_Format$Type == "FLOAT") {
    ImgSubset <- readBin(fid, numeric(), n = as.double(nbRows) * as.double(HDR$samples) * as.double(HDR$bands), size = Image_Format$Bytes, signed = Image_Format$Signed, endian = Image_Format$ByteOrder)
  }
  close(fid)
  # reshape ImgSubset in a 2D matrix in order to get selected pixels based on index
  ImgSubset <- array(aperm(array(ImgSubset, dim = c(HDR$samples, HDR$bands, nbRows)), c(3, 1, 2)), c(HDR$samples * nbRows, HDR$bands))
  # coordinates of the samples correspond to the whole image: need to convert to fit sample
  coordPix_List[, 1] <- coordPix_List[, 1] - minRow + 1
  # Get index of each pixel
  IndPix <- (coordPix_List[, 2] - 1) * nbRows + coordPix_List[, 1]
  # Get samples from subset
  Sample_Sel <- ImgSubset[IndPix, ]
  rm(ImgSubset)
  gc()
  return(Sample_Sel)
}

#' extracts pixels from image based on their coordinates
#'
#' @param ImPath character. path for image
#' @param coordPix numeric. pixel coordinates
#' @param MaxRAM numeric.
#' @param Already.Split Boolean.
#'
#' @return samples from image corresponding to coordPix
#' @export

extract_samples_from_image <- function(ImPath, coordPix, MaxRAM = FALSE, Already.Split = FALSE) {
  # get image header data
  ImPathHDR <- get_HDR_name(ImPath)
  HDR <- read_ENVI_header(ImPathHDR)

  # compute the ranking of initial pixel list compared to index ranking
  if (typeof(coordPix) == "double" ){
    if (dim(coordPix)[2] == 2) {
      if (dim(coordPix)[1] >= 2) {
        coordPix_tmp <- list()
        coordPix_tmp$row <- coordPix[, 1]
        coordPix_tmp$col <- coordPix[, 2]
      } else if (dim(coordPix)[1] == 1) {
        coordPix_tmp <- list()
        coordPix_tmp$row <- coordPix[1]
        coordPix_tmp$col <- coordPix[2]
      }
    }
  } else if (typeof(coordPix) == "list" & length(grep("row", names(coordPix))) > 0 & length(grep("col", names(coordPix))) > 0) {
    coordPix_tmp <- coordPix
  }
  # initial index value of the pixels requested in the image, following original ranking
  Index_Init <- sub2ind(HDR, coordPix_tmp)
  # rank of the initial list of pixels
  Rank_Index_Init <- sort(Index_Init, index = TRUE)
  # convert
  if (typeof(coordPix) == "list" & length(grep("row", names(coordPix))) > 0 & length(grep("col", names(coordPix))) > 0) {
    coordPix <- cbind(coordPix$row, coordPix$col)
  }
  # divide data access based on the size of the image (to avoid reading 30 Gb file at once...)
  Image_Format <- ENVI_type2bytes(HDR)
  lenTot <- as.double(HDR$lines) * as.double(HDR$samples) * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image_Format$Bytes) / (1024^3)
  # maximum image size read at once. If image larger, then reads in multiple pieces
  if (MaxRAM == FALSE) {
    LimitSizeGb <- 0.25
  } else {
    LimitSizeGb <- MaxRAM
  }
  if (ImSizeGb < LimitSizeGb | Already.Split == TRUE) {
    Lines_Per_Read <- HDR$lines
    nbPieces <- 1
  } else {
    # nb of lines corresponding to LimitSizeGb
    OneLine <- as.double(HDR$samples) * as.double(HDR$bands) * Image_Format$Bytes
    Lines_Per_Read <- floor(LimitSizeGb * (1024^3) / OneLine)
    # number of pieces to split the image into
    nbPieces <- ceiling(HDR$lines / Lines_Per_Read)
  }

  # here split based on even number of pxiels to sample
  # should be based on image size in order to avoid memory
  # problems if target pixels  unevenly distributed in image
  coordPix_List <- split_pixel_samples(coordPix, Lines_Per_Read)
  Sample_Sel <- list()
  for (i in 1:length(coordPix_List)) {
    Sample_Sel[[i]] <- extract_pixels(coordPix_List[[i]], ImPath, HDR)
  }
  Sample_Sel <- do.call("rbind", Sample_Sel)
  coordPix_List <- do.call("rbind", coordPix_List)

  # re-sort samples
  Sample_Sort <- list()
  Sample_Sort$row <- coordPix_List[, 1]
  Sample_Sort$col <- coordPix_List[, 2]
  Coord_Pixels <- sub2ind(HDR, Sample_Sort)
  # rank of the pixels extracted from the image
  Rank_Pixels <- sort(Coord_Pixels, index = TRUE)
  # pixel re-ordering needs to be performed in order to get back to Rank_Index_Init
  # first re-order pixels in order to follow increasing index value
  # Sample_Sort$index   = Coord_Pixels[Rank_Pixels$ix]
  # Sample_Sort$row     = Sample_Sort$row[Rank_Pixels$ix]
  # Sample_Sort$col  = Sample_Sort$col[Rank_Pixels$ix]
  # if bug, check this line
  # Sample_Sel          = Sample_Sel[Rank_Pixels$ix]
  Sample_Sel <- Sample_Sel[Rank_Pixels$ix, ]

  # then apply initial ranking as defined by Rank_Index_Init
  if (dim(Sample_Sel)[1] > 1) {
    Sample_Sel[Rank_Index_Init$ix, ] <- Sample_Sel
  } else if (dim(Sample_Sel)[1] == 1) {
    Sample_Sel <- Sample_Sel
  }
  return(Sample_Sel)
}

#' Extract bands of sparse pixels in image data cube
#' @param ImPath character. Path to the image cube
#' @param rowcol matrix or data.frame with two columns: row, col.
#' If columns are not named, 1st=row, 2nd=col.
#' @param MaxRAM numeric. Maximum memory use at block reading.
#' It constrains the maximum number rows of a block
#'
#' @return matrix. Rows are corresponding to the samples, columns are the bands.
#' @importFrom raster brick
#' @import stars
#' @export

extract.big_raster <- function(ImPath, rowcol, MaxRAM=.50){

  if(!is.data.frame(rowcol)){
    rowcol <- as.data.frame(rowcol)
  }

  if(!all(c('row', 'col') %in% colnames(rowcol))){
    warning('Columns row,col not found in rowcol argument. The two first columns are considered as row, col respectively.')
    colnames(rowcol)[1:2]= c('row', 'col')
  }

  metarast <- raster(ImPath)
  # updated raster package: do not use brick with 2D raster
  if (nbands(metarast)>1){
    rasterInfo <- raster::brick(ImPath)
  } else{
    rasterInfo <- metarast
  }

  # nbytes = as.numeric(substring(dataType(rasterInfo), 4, 4))
  # stars converts automatically values to numeric
  nbytes <- 8
  ImgSizeGb <- prod(dim(rasterInfo))*nbytes/2^30
  LineSizeGb <- prod(dim(rasterInfo)[2:3])*nbytes/2^30
  LinesBlock <- floor(MaxRAM/LineSizeGb)
  rowcol$rowInBlock <- ((rowcol$row-1) %% LinesBlock)+1  # row in block number
  rowcol$block <- floor((rowcol$row-1)/LinesBlock)+1  # block number
  rowcol$sampleIndex <- 1:nrow(rowcol)  # sample index to reorder result

  sampleList = lapply(unique(rowcol$block), function(iblock){
    rc <- rowcol[rowcol$block==iblock,]
    rr <- range(rc$row)
    nYSize <- diff(rr)+1
    nXSize <- max(rc$col)
    # stars data cube dimension order is x*y*band
    ipix_stars <- (rc$rowInBlock-min(rc$rowInBlock))*nXSize+rc$col
    # get driver
    driver <- attr(rgdal::GDALinfo(ImPath,returnStats = FALSE), 'driver')
    values <- read_stars(ImPath, RasterIO =list(nXSize=nXSize, nYOff=rr[1], nYSize=nYSize),proxy = FALSE, driver=driver)[[1]]
    values <- matrix(values, nrow=nYSize*nXSize)
    res <- cbind(rc$sampleIndex, values[ipix_stars, ])
    rm('values')
    gc()
    return(res)
  })

  samples = do.call(rbind, sampleList)
  samples = samples[order(samples[,1]),2:ncol(samples)]

  return(samples)
}

#' extract random subset of pixels from an image
#'
#' @param ImPath character. path for the image to sample
#' @param MaskPath character. path for the corresponding mask
#' @param nb_partitions numeric. number of repetitions of kmeans
#' @param Pix_Per_Partition numeric.
#' @param kernel numeric.
#' @param MaxRAM numeric.
#'
#' @importFrom matlab ones
#' @import raster
#' @importFrom mmand erode
#' @importFrom data.table data.table rbindlist setorder
#' @importFrom matrixStats rowAnys
#' @return list, see Details
#' @details
#' The returned list contains:
#' - DataSubset: matrix of NxP of N samples and P bands
#' - nbPix2Sample: integer giving the number of pixels sampled (only central pixel of kernel)
#' - coordPix: a data.table with columns 'row', 'col' of pixel in the image corresponding to each row of DataSubset, and if kernel is not NULL
#' Kind (Kernel index) and 'id' the sample ID to be used with the kernel
#' @export

get_random_subset_from_image <- function(ImPath, MaskPath, nb_partitions, Pix_Per_Partition, kernel=NULL,MaxRAM = 0.5) {

  metarast <- raster(ImPath)
  # updated raster package: do not use brick with 2D raster
  if (nbands(metarast)>1){
    rasterInfo <- raster::brick(ImPath)
  } else{
    rasterInfo <- metarast
  }
  nbPix2Sample <- nb_partitions * Pix_Per_Partition
  # get total number of pixels
  rdim <- dim(rasterInfo)
  nlines <- rdim[1]
  nsamples <- rdim[2]
  nbpix <- ncell(rasterInfo)
  # 1- Exclude masked pixels from random subset
  # Read Mask
  if ((!MaskPath == "") & (!MaskPath == FALSE)) {
    mask <- matrix(t(raster(MaskPath)),ncol= nsamples,nrow = nlines)
  } else {
    mask <- array(1, dim = c(nlines, nsamples))
  }

  if(is.matrix(kernel)){
    # erode mask with kernel, to keep valid central pixels and neighbours
    mask = matlab::padarray(mask, c(1,1), padval=0, direction='both')
    mask = erode(mask, (kernel!=0)*1)
    mask = mask[2:(nrow(mask)-1), 2:(ncol(mask)-1)]
  }

  # get a list of samples among unmasked pixels

  ValidPixels <- which(mask > 0)
  NbValidPixels <- length(ValidPixels)
  # Check if number of valid pixels is compatible with number of pixels to be extracted
  # if number of pixels to sample superior to number of valid pixels, then adjust iterations
  if (NbValidPixels < nbPix2Sample) {
    nbPix2Sample <- NbValidPixels
    nb_partitions <- ceiling(NbValidPixels / Pix_Per_Partition)
    Pix_Per_Partition <- floor(NbValidPixels / nb_partitions)
    nbPix2Sample <- nb_partitions * Pix_Per_Partition
  }

  # Select a random subset of nbPix2Sample
  seed <- sample(10000)[1]
  set.seed(seed)
  pixselected <- sample(ValidPixels, nbPix2Sample)
  Row <- ((pixselected - 1) %% nlines) + 1 # row
  Column <- floor((pixselected - 1) / nlines) + 1 # column
  if(is.matrix(kernel)){
    coordPixK = list()
    mesh=matlab::meshgrid(-(ncol(kernel)%/%2):(ncol(kernel)%/%2), -(nrow(kernel)%/%2):(nrow(kernel)%/%2))
    for(p in which(kernel!=0)){
      coordPixK[[p]] = data.table(row = Row+mesh$y[p], col = Column+mesh$x[p], id=1:length(Row))
    }
    coordPix = rbindlist(coordPixK, idcol='Kind')
    # Order along coordPix$id for further use in noise, mnf
    setorder(coordPix, 'id')
  }else{
    coordPix = data.table(row = Row, col = Column, id = 1:length(Row))
  }
  # sort based on .bil dim order, i.e. band.x.y or band.col.row
  # TODO: sorting may not be necessary anymore, neither unique coordinates
  # setorder(coordPix, col, row)

  # make unique
  ucoordPix <- unique(coordPix[,c('row','col')])
  ucoordPix[['sampleIndex']] = 1:nrow(ucoordPix)

  # 2- Extract samples from image
  # TODO: if coordPix can be a data.frame it would be easier
  # Sample_Sel <- biodivMapR:::extract_samples_from_image(ImPath, ucoordPix)
  Sample_Sel <- extract.big_raster(ImPath, ucoordPix[,1:2])
  samplePixIndex <- merge(coordPix, ucoordPix, by=c('row', 'col'), sort = FALSE)
  # randomize! It should already be random from the sample operation on mask valid pixels
  Sample_Sel <- Sample_Sel[samplePixIndex$sampleIndex, ]
  samplePixIndex[['sampleIndex']]=NULL

  # remove NA and Inf pixels
  if(any(is.na(Sample_Sel) | is.infinite(Sample_Sel))){
    print('Removing pixels with NA values.')
    rmrows <- rowAnys(is.na(Sample_Sel) | is.infinite(Sample_Sel))
    rmpix <- unique(samplePixIndex$id[rmrows])
    Sample_Sel = Sample_Sel[!(samplePixIndex$id %in% rmpix),]
    samplePixIndex = samplePixIndex[!(samplePixIndex$id %in% rmpix)]
    nbPix2Sample = length(unique(samplePixIndex$id))
  }

  ### FLORIAN: REMOVED RANDOMIZATION AS IT SEEMED USELESS
  # Sample_Sel <- Sample_Sel[sample(dim(Sample_Sel)[1]), ]
  # coordPix <- coordPix[sample(dim(Sample_Sel)[1]), ]
  # TODO: check how returned coordPix is used, as it is now a data.frame
  my_list <- list("DataSubset" = Sample_Sel, "nbPix2Sample" = nbPix2Sample,"coordPix"=samplePixIndex)
  return(my_list)
}

# does the system work with little endians or big endians?
#
# @return ByteOrder
#' @import tools
get_byte_order <- function() {
  if (.Platform$endian == "little") {
    ByteOrder <- 0
  } else if (.Platform$endian == "big") {
    ByteOrder <- 1
  }
  return(ByteOrder)
}

#' get hdr name from image file name, assuming it is BIL format
#'
#' @param ImPath character. ath of the image
#' @param showWarnings boolean. set TRUE if warning because HDR does not exist
#'
#' @return corresponding hdr
#' @import tools
#' @export
get_HDR_name <- function(ImPath,showWarnings=TRUE) {
  if (file_ext(ImPath) == "") {
    ImPathHDR <- paste(ImPath, ".hdr", sep = "")
  } else if (file_ext(ImPath) == "bil") {
    ImPathHDR <- gsub(".bil", ".hdr", ImPath)
  } else if (file_ext(ImPath) == "zip") {
    ImPathHDR <- gsub(".zip", ".hdr", ImPath)
  } else {
    ImPathHDR <- paste(file_path_sans_ext(ImPath), ".hdr", sep = "")
  }
  if (showWarnings==TRUE){
    if (!file.exists(ImPathHDR)) {
      message("WARNING : COULD NOT FIND HDR FILE")
      print(ImPathHDR)
      message("Process may stop")
    }
  }
  return(ImPathHDR)
}

# gets rank of spectral bands in an image
#
# @param Spectral_Bands wavelength (nm) of the spectral bands to be found
# @param wavelength wavelength (nm) of all wavelengths in the image
#
# @return rank of all spectral bands of interest in the image and corresponding wavelength

get_image_bands <- function(Spectral_Bands, wavelength) {
  ImBand <- c()
  Distance2WL <- c()
  for (band in Spectral_Bands) {
    Closest_Band <- order(abs(wavelength - band))[1]
    ImBand <- c(ImBand, Closest_Band)
    Distance2WL <- c(Distance2WL, abs(wavelength[Closest_Band] - band))
  }
  my_list <- list("ImBand" = ImBand, "Distance2WL" = Distance2WL)
  return(my_list)
}

#' extract coordinates of the footprint of elements of a shapefile in a raster
#'
#' @param path_SHP character. path for the shapefile
#' @param path_Raster character. path for the raster
#' @param IDshp character. field in the shapefile determining ID of element
#'
#' @return list. vector_coordinates and vector_ID for each element in the vector file
#' @importFrom tools file_path_sans_ext
#' @importFrom rgdal readOGR
#' @import raster
#' @export

get_polygonCoord_from_Shp <- function(path_SHP,path_Raster,IDshp=NULL){
  # prepare for possible reprojection
  Dir_Vector <- dirname(path_SHP)
  Name_Vector <- tools::file_path_sans_ext(basename(path_SHP))
  print(paste('Reading pixels coordinates for polygons in ',Name_Vector,sep=''))
  # File.Vector.reproject <- paste(Dir_Vector.reproject,'/',Name_Vector,'.shp','sep'='')
  if (file.exists(paste(file_path_sans_ext(path_SHP),'.shp',sep=''))){
    Plot <- rgdal::readOGR(Dir_Vector,Name_Vector,verbose = FALSE)
    # check if vector and rasters are in the same referential
    # if not, convert vector file
    if (!raster::compareCRS(raster(path_Raster), Plot)){
      stop('Raster and Plots have different projection. Plots should be reprojected to Raster CRS')
    }
  } else if (file.exists(paste(path_SHP,'kml','sep'='.'))){
    print('Please convert vector file to shpfile')
  }
  # extract data corresponding to the Raster_SpectralSpecies
  vector_coordinates <- extract_pixels_coordinates.From.OGR(path_Raster,Plot)
  # if the different elements can be identified
  vector_ID <- c()
  for (elem in 1:length(vector_coordinates)){
    if (!is.null(IDshp)){
      vector_ID <- c(vector_ID,Plot[[IDshp]])
    } else {
      vector_ID <- c(vector_ID,'No ID Available')
    }
  }
  my_list <- list("vector_coordinates" = vector_coordinates, "vector_ID" = vector_ID)
  return(my_list)
}

#' convert image coordinates from index to X-Y
#'
#' @param Raster image raster object
#' @param Image_Index coordinates corresponding to the raster
#' @export

ind2sub <- function(Raster, Image_Index) {
  c <- ((Image_Index - 1) %% Raster@ncols) + 1
  r <- floor((Image_Index - 1) / Raster@ncols) + 1
  my_list <- list("col" = c, "row" = r)
  return(my_list)
}

# convert image coordinates from index to X-Y
# image coordinates are given as index = (ID.col-1) * total.lines + ID.row
#
# @param Raster image raster object
# @param Image_Index coordinates corresponding to the raster
ind2sub2 <- function(Raster, Image_Index) {
  r <- ((Image_Index - 1) %% Raster@nrows) + 1
  c <- floor((Image_Index - 1) / Raster@nrows) + 1
  my_list <- list("col" = c, "row" = r)
  return(my_list)
}

#' This function computes interquartile range (IQR) criterion, which can be used
#' as a criterion for outlier detection
#'
#' @param DistVal numeric. vector of distribution of values
#' @param weightIRQ numeric. weighting factor appplied to IRQ to define lower and upper boudaries for outliers
#'
#' @return outlier_IQR numeric. band numbers of original sensor corresponding to S2
#' @importFrom stats IQR quantile
#' @export
IQR_outliers <- function(DistVal,weightIRQ = 1.5){
  iqr <- IQR(DistVal, na.rm=TRUE)
  range_IQR <- c(quantile(DistVal, 0.25,na.rm=TRUE),quantile(DistVal, 0.75,na.rm=TRUE))
  outlier_IQR <- c(range_IQR[1]-weightIRQ*iqr,range_IQR[2]+weightIRQ*iqr)
  return(outlier_IQR)
}

#' applies mean filter to an image
#'
#' @param Image numeric. image defined as 2d matrix
#' @param SizeFilt numeric. size of the window to compute mean filter
#' @param NA_remove boolean. Should NA be removed?
#'
#' @return rank of all spectral bands of interest in the image and corresponding wavelength
#' @importFrom matlab padarray
#' @export

mean_filter <- function(Image, SizeFilt,NA_remove = FALSE) {
  nbi <- dim(Image)[1]
  nbj <- dim(Image)[2]
  # if matrix: convert into array
  if (is.na(dim(Image)[3])){
    Image <- array(Image,dim = c(nbi,nbj,1))
  }
  ImageSmooth <- array(0,c(nbi,nbj,dim(Image)[3]))
  for (band in 1: dim(Image)[3]){
    E <- padarray(Image[,,band], c(SizeFilt, SizeFilt), "symmetric", "both")
    ImageSmooth_tmp <- matrix(0, nrow = nbi, ncol = nbj)
    Mat2D <- MatSun <- matrix(0, nrow = ((2 * SizeFilt) + 1)^2, ncol = nbj)
    spl <- split(1:nbj, 1:((2 * SizeFilt) + 1))
    mid <- ceiling((((2 * SizeFilt) + 1)^2) / 2)
    for (i in (SizeFilt + 1):(nbi + SizeFilt)) {
      for (j in 1:((2 * SizeFilt) + 1)) {
        # create a 2D matrix
        Mat2D[, spl[[j]]] <- matrix(E[(i - SizeFilt):(i + SizeFilt), (spl[[j]][1]):(tail(spl[[j]], n = 1) + 2 * SizeFilt)], nrow = ((2 * SizeFilt) + 1)^2)
      }
      ImageSmooth_tmp[(i - SizeFilt), ] <- colMeans(Mat2D, na.rm = TRUE)
    }
    if (NA_remove == TRUE){
      ImageSmooth_tmp[which(is.na(Image[,,band]))] <- NA
    }
    ImageSmooth[,,band] <- ImageSmooth_tmp
  }
  return(ImageSmooth)
}

#' reads a subset from a binary image
#'
#' @param Byte_Start numeric. location of byte where to start reading in the image
#' @param nbLines numeric. number of lines to read
#' @param lenBin numeric. number of elements to read
#' @param ImPath character. path for the image
#' @param ImBand numeric. bands of interest
#' @param jpix numeric. number of columns in the image
#' @param nbChannels numeric. total number of channels in the image
#' @param Image_Format character. type of data (INT/FLOAT)
#'
#' @return data corresponding to the subset in original 3D format
#' @export

read_bin_subset <- function(Byte_Start, nbLines, lenBin, ImPath, ImBand, jpix, nbChannels, Image_Format) {
  # number of bands to be kept
  nbSubset <- length(ImBand)
  # open file
  fid <- file(
    description = ImPath, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  if (!Byte_Start == 1) {
    # skip the beginning of the file of not wanted
    seek(fid, where = as.double(Image_Format$Bytes * (Byte_Start - 1)), origin = "start", rw = "read")
  }
  if (Image_Format$Type == "INT") {
    linetmp <- readBin(fid, integer(), n = lenBin, size = Image_Format$Bytes, endian = Image_Format$ByteOrder)
  } else if (Image_Format$Type == "FLOAT") {
    linetmp <- readBin(fid, numeric(), n = lenBin, size = Image_Format$Bytes, endian = Image_Format$ByteOrder)
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

#' Reads ENVI hdr file
#'
#' @param HDRpath Path of the hdr file
#'
#' @return list of the content of the hdr file
#' @export

read_ENVI_header <- function(HDRpath) {
  # header <- paste(header, collapse = "\n")
  if (!grepl(".hdr$", HDRpath)) {
    stop("File extension should be .hdr")
  }
  HDR <- readLines(HDRpath)
  ## check ENVI at beginning of file
  if (!grepl("ENVI", HDR[1])) {
    stop("Not an ENVI header (ENVI keyword missing)")
  } else {
    HDR <- HDR [-1]
  }
  ## remove curly braces and put multi-line key-value-pairs into one line
  HDR <- gsub("\\{([^}]*)\\}", "\\1", HDR)
  l <- grep("\\{", HDR)
  r <- grep("\\}", HDR)

  if (length(l) != length(r)) {
    stop("Error matching curly braces in header (differing numbers).")
  }

  if (any(r <= l)) {
    stop("Mismatch of curly braces in header.")
  }

  HDR[l] <- sub("\\{", "", HDR[l])
  HDR[r] <- sub("\\}", "", HDR[r])

  for (i in rev(seq_along(l))) {
    HDR <- c(
      HDR [seq_len(l [i] - 1)],
      paste(HDR [l [i]:r [i]], collapse = "\n"),
      HDR [-seq_len(r [i])]
    )
  }

  ## split key = value constructs into list with keys as names
  HDR <- sapply(HDR, split_line, "=", USE.NAMES = FALSE)
  names(HDR) <- tolower(names(HDR))

  ## process numeric values
  tmp <- names(HDR) %in% c(
    "samples", "lines", "bands", "header offset", "data type",
    "byte order", "default bands", "data ignore value",
    "wavelength", "fwhm", "data gain values"
  )
  HDR [tmp] <- lapply(HDR [tmp], function(x) {
    as.numeric(unlist(strsplit(x, ",")))
  })

  return(HDR)
}

#' read specific image bands from image
#'
#' @param ImPath Path of the image to read
#' @param HDR Header for the image
#' @param ImBand Bands to be read
#'
#' @return Image_Subset information corresponding to ImBand
#' @import stars

read_image_bands <- function(ImPath, HDR, ImBand) {
  # first get image format
  Image_Format <- ENVI_type2bytes(HDR)
  ipix <- HDR$lines
  jpix <- HDR$samples
  nbChannels <- HDR$bands
  nbSubset <- length(ImBand)
  Image_Subset <- array(NA,c(ipix,jpix,length(ImBand)))
  i <- 0
  for (band in ImBand){
    i <- i+1
    bndtmp = t(matrix(raster(ImPath, band = ImBand[i]),nrow = HDR$samples,ncol = HDR$lines))
    Image_Subset[,,i] <- array(bndtmp,c(ipix,jpix,1))
  }
  return(Image_Subset)
}

#' reads subset of lines from an image
#'
#' @param ImPath path for the image
#' @param HDR header information corresponding to the image
#' @param Line_Start which line to start reading
#' @param Lines_To_Read number of lines to read
#' @param ImgFormat character. is it a matrix (2D) or a raster (3D)?
#'
#' @return data corresponding to the subset in original 3D format
#' @export

read_image_subset <- function(ImPath, HDR, Line_Start,Lines_To_Read,ImgFormat='3D'){
  # list of pixels to be extracted
  ListRows <- seq(Line_Start,Line_Start+Lines_To_Read-1)
  ListCols <- seq(1,HDR$samples)
  # list pixels
  coordPix <- expand.grid(ListRows,ListCols)
  names(coordPix) <- c('row','col')
  coordPix[['sampleIndex']] = 1:nrow(coordPix)
  # 2- Extract samples from image
  Sample_Sel <- extract.big_raster(ImPath, coordPix[,1:2])
  Sample_Sel <- array(matrix(as.numeric(unlist(Sample_Sel)),ncol = HDR$bands),c(Lines_To_Read,HDR$samples,HDR$bands))
  if (ImgFormat == "2D") {
    Sample_Sel <- array(Sample_Sel, c(Lines_To_Read * HDR$samples, HDR$bands))
  }
  if (ImgFormat == "Shade") {
    Sample_Sel <- matrix(Sample_Sel, Lines_To_Read, HDR$samples, byrow = F)
  }
  return(Sample_Sel)
}

#' reads subset of an ENVI BIL image
#'
#' @param ImPath path for the image
#' @param HDR header information corresponding to the image
#' @param Byte_Start location of byte where to start reading in the image
#' @param lenBin number of elements to read
#' @param nbLines number of lines to read
#' @param Image_Format type of data (INT/FLOAT)
#' @param ImgFormat should output be 2D or 3D (original image format)?
#'
#' @return data corresponding to the subset in original 3D format
#' @export

read_BIL_image_subset <- function(ImPath, HDR, Byte_Start, lenBin, nbLines, Image_Format, ImgFormat) {
  fidIm <- file(
    description = ImPath, open = "rb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )

  # read relevant portion
  if (!Byte_Start == 1) {
    # skip the beginning of the file of not wanted
    seek(fidIm, where = Image_Format$Bytes * (Byte_Start - 1), origin = "start", rw = "read")
  }
  if (Image_Format$Type == "INT") {
    if (HDR$`data type` == 1) {
      Image_Chunk <- readBin(fidIm, integer(), n = lenBin, size = Image_Format$Bytes, signed = FALSE, endian = Image_Format$ByteOrder)
    } else {
      Image_Chunk <- readBin(fidIm, integer(), n = lenBin, size = Image_Format$Bytes, signed = TRUE, endian = Image_Format$ByteOrder)
    }
  } else if (Image_Format$Type == "FLOAT") {
    Image_Chunk <- readBin(fidIm, numeric(), n = lenBin, size = Image_Format$Bytes, endian = Image_Format$ByteOrder)
  }

  if (ImgFormat == "3D" | ImgFormat == "2D") {
    Image_Chunk <- aperm(array(Image_Chunk, dim = c(HDR$samples, HDR$bands, nbLines)), c(3, 1, 2))
  }
  if (ImgFormat == "2D") {
    Image_Chunk <- array(Image_Chunk, c(nbLines * HDR$samples, HDR$bands))
  }
  if (ImgFormat == "Shade") {
    Image_Chunk <- matrix(Image_Chunk, nbLines, HDR$samples, byrow = T)
  }
  close(fidIm)
  return(Image_Chunk)
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

split_line <- function(x, separator, trim.blank = TRUE) {
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

#' splits a set of pixels to be sampled in an image based on number of lines, not number of samples
#'
#' @param coordPix numeric. coordinates of pixels to sample
#' @param Lines_Per_Read numeric. max number of lines per read for memory concerns
#'
#' @return coordPix_List list of pixel coordinates
#' @export

split_pixel_samples <- function(coordPix, Lines_Per_Read) {
  # maximum Lines_Per_Read
  if (dim(coordPix)[1] > 1) {
    nb.Lines <- max(coordPix[, 1]) - min(coordPix[, 1]) + 1
  } else if (dim(coordPix)[1] == 1) {
    nb.Lines <- max(coordPix[1]) - min(coordPix[1]) + 1
  }
  nb.Pieces <- ceiling(nb.Lines / Lines_Per_Read)

  if (dim(coordPix)[1] > 1) {
    Min.Line <- ceiling(seq(min(coordPix[, 1]), max(coordPix[, 1]), length.out = nb.Pieces + 1))
  } else if (dim(coordPix)[1] == 1) {
    Min.Line <- ceiling(seq(min(coordPix[1]), max(coordPix[1]), length.out = nb.Pieces + 1))
  }

  Max.Line <- Min.Line - 1
  Min.Line <- Min.Line[-length(Min.Line)]
  Max.Line <- Max.Line[-1]
  Max.Line[length(Max.Line)] <- Max.Line[length(Max.Line)] + 1
  coordPix_List <- list()
  ii <- 0
  for (i in 1:length(Min.Line)) {
    if (dim(coordPix)[1] > 1) {
      selpix <- which(coordPix[, 1] >= Min.Line[i] & coordPix[, 1] <= Max.Line[i])
    } else if (dim(coordPix)[1] == 1) {
      selpix <- which(coordPix[1] >= Min.Line[i] & coordPix[1] <= Max.Line[i])
    }
    if (length(selpix) > 1) {
      ii <- ii + 1
      coordPix_List[[ii]] <- coordPix[selpix, ]
    } else if (length(selpix) == 1) {
      ii <- ii + 1
      coordPix_List[[ii]] <- coordPix
    }
  }
  return(coordPix_List)
}

#' updates an existing mask
#'
#' @param MaskPath character. path for original mask (may not exist)
#' @param HDR list. header correpondingproviding general info about data format
#' @param Mask numeric. data to be used in the mask
#' @param MaskPath_Update character. path for the updated mask
#'
#' @return MaskPath_Update
#' @export

update_shademask <- function(MaskPath, HDR, Mask, MaskPath_Update) {
  ipix <- HDR$lines
  jpix <- HDR$samples
  nbpix <- as.double(ipix) * as.double(jpix)
  # if MaskPath provided
  if ((!MaskPath == "") & (!MaskPath == FALSE)) {
    # read shade mask
    fid <- file(
      description = MaskPath, open = "rb", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    lenBin <- nbpix
    MaskTmp <- readBin(fid, integer(), n = lenBin, size = 1)
    close(fid)
    MaskTmp <- aperm(array(MaskTmp, dim = c(jpix, ipix)))
    # multiply by Mask
    Mask <- MaskTmp * Mask
  }
  Mask <- array(Mask, c(ipix, jpix, 1))
  Mask <- aperm(Mask, c(2, 3, 1))
  fidOUT <- file(
    description = MaskPath_Update, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  writeBin(c(as.integer(Mask)), fidOUT, size = 1, endian = .Platform$endian, useBytes = FALSE)
  close(fidOUT)
  # write updated shademask
  HDR_Update <- HDR
  HDR_Update$description <- "Mask produced from radiometric filtering"
  HDR_Update$bands <- 1
  HDR_Update$`data type` <- 1
  HDR_Update$`file type` <- NULL
  HDR_Update$`band names` <- "Mask"
  HDR_Update$`default stretch` <- '0 1 linear'
  HDR_Update$wavelength <- NULL
  HDR_Update$fwhm <- NULL
  HDR_Update$resolution <- NULL
  HDR_Update$bandwidth <- NULL
  HDR_Update$purpose <- NULL
  HDR_Update$`default bands` <- NULL
  HDR_Update$`data gain values` <- NULL
  HDRpath <- paste(MaskPath_Update, ".hdr", sep = "")
  write_ENVI_header(HDR_Update, HDRpath)
  return(MaskPath_Update)
}

#' defines which byte should be read for each part of an image split in nbPieces
#'
#' @param HDR list. header info
#' @param nbPieces numeric. number of pieces resulting from image split
#'
#' @return location of the bytes corresponding to beginning and end of each piece, and corresponding number of lines
#' @export

where_to_read <- function(HDR, nbPieces) {
  Data_Per_Line <- as.double(HDR$samples) * as.double(HDR$bands)
  lenTot <- as.double(HDR$lines) * Data_Per_Line
  # starting line for each chunk
  Start_Per_Chunk <- ceiling(seq(1, (HDR$lines + 1), length.out = nbPieces + 1))
  # elements in input data
  lb <- 1 + ((Start_Per_Chunk - 1) * Data_Per_Line)
  ub <- lb - 1
  ReadByte_Start <- lb[1:nbPieces]
  ReadByte_End <- ub[2:(nbPieces + 1)]
  Lines_Per_Chunk <- diff(Start_Per_Chunk)
  my_list <- list("ReadByte_Start" = ReadByte_Start, "ReadByte_End" = ReadByte_End, "Lines_Per_Chunk" = Lines_Per_Chunk,'Line_Start'=Start_Per_Chunk)
  return(my_list)
}

#' defines which byte should be read for each part of an image split in nbPieces
#'
#' @param HDR list. header info
#' @param nbPieces numeric. number of pieces resulting from image split
#' @param SE_Size numeric. size of structuring element (window)
#'
#' @return location of the bytes corresponding to beginning and end of each piece, and corresponding number of lines
#' @export

where_to_read_kernel <- function(HDR, nbPieces, SE_Size) {
  Data_Per_Line <- as.double(HDR$samples) * as.double(HDR$bands)
  lenTot <- as.double(HDR$lines) * Data_Per_Line
  # starting line for each chunk
  Start_Per_Chunk <- ceiling(seq(1, (HDR$lines + 1), length.out = nbPieces + 1))
  Start_Per_Chunk <- Start_Per_Chunk - Start_Per_Chunk %% SE_Size + 1
  # elements in input data
  lb <- 1 + ((Start_Per_Chunk - 1) * Data_Per_Line)
  ub <- lb - 1
  ReadByte_Start <- lb[1:nbPieces]
  ReadByte_End <- ub[2:(nbPieces + 1)]
  Lines_Per_Chunk <- diff(Start_Per_Chunk)
  my_list <- list("ReadByte_Start" = ReadByte_Start, "ReadByte_End" = ReadByte_End, "Lines_Per_Chunk" = Lines_Per_Chunk)
  return(my_list)
}

#' defines which byte should be written for each part of an image split in nbPieces
#'
#' @param HDR_SS list. header info for SpectralSpecies file
#' @param HDR_SSD list. header info for SpectralSpecies_Distribution file
#' @param nbPieces numeric. number of pieces resulting from image split
#' @param SE_Size numeric.
#'
#' @return location of the bytes corresponding to beginning and end of each piece, and corresponding number of lines
#' @export

where_to_write_kernel <- function(HDR_SS, HDR_SSD, nbPieces, SE_Size) {
  Data_Per_Line_SS <- as.double(HDR_SS$samples) * as.double(HDR_SS$bands)
  Data_Per_Line_SSD <- as.double(HDR_SSD$samples) * as.double(HDR_SSD$bands)

  # starting line for each chunk of spectral species
  Start_Per_Chunk <- ceiling(seq(1, (HDR_SS$lines + 1), length.out = nbPieces + 1))
  Start_Per_Chunk <- Start_Per_Chunk - Start_Per_Chunk %% SE_Size
  Start_Per_Chunk.SSD <- (Start_Per_Chunk / SE_Size) + 1

  # elements in input data
  Image_Format <- ENVI_type2bytes(HDR_SSD)
  lb_SSD <- 1 + (((Start_Per_Chunk.SSD - 1) * Data_Per_Line_SSD) * Image_Format$Bytes)
  ub_SSD <- lb_SSD - 1
  ReadByte_Start.SSD <- lb_SSD[1:nbPieces]
  ReadByte_End.SSD <- ub_SSD[2:(nbPieces + 1)]
  Lines_Per_Chunk.SSD <- diff(Start_Per_Chunk.SSD)

  my_list <- list("ReadByte_Start" = ReadByte_Start.SSD, "ReadByte_End" = ReadByte_End.SSD, "Lines_Per_Chunk" = Lines_Per_Chunk.SSD)
  return(my_list)
}

#' writes ENVI hdr file
#'
#' @param HDR content to be written
#' @param HDRpath Path of the hdr file
#'
#' @return
#' @importFrom stringr str_count
#' @export

write_ENVI_header <- function(HDR, HDRpath) {
  h <- lapply(HDR, function(x) {
    if (length(x) > 1 || (is.character(x) && str_count(x, "\\w+") > 1)) {
      x <- paste0("{", paste(x, collapse = ","), "}")
    }
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c("ENVI", paste(names(HDR), h, sep = " = ")), con = HDRpath)
  return("")
}

#' write an image which size is > 2**31-1
#'
#' @param ImgWrite numeric. Image as array
#' @param ImagePath character. Path where the image should be written
#' @param HDR list. Image header
#' @param Image_Format list. description of data format corresponding to ENVI type
#'
#' @return None
#' @export

Write_Big_Image <- function(ImgWrite,ImagePath,HDR,Image_Format){
  nbPieces <- split_image(HDR, LimitSizeGb = 1.8)
  SeqRead_Image <- where_to_read(HDR, nbPieces)
  fidOUT <- file(
    description = ImagePath, open = "wb", blocking = TRUE,
    encoding = getOption("encoding"), raw = FALSE
  )
  close(fidOUT)
  # for each piece of image
  for (i in 1:nbPieces) {
    print(paste("Writing Image, piece #", i, "/", nbPieces))
    # read image and mask data
    Byte_Start <- SeqRead_Image$ReadByte_Start[i]
    Line_Start <- SeqRead_Image$Line_Start[i]
    nbLines <- SeqRead_Image$Lines_Per_Chunk[i]
    lenBin <- SeqRead_Image$ReadByte_End[i] - SeqRead_Image$ReadByte_Start[i] + 1
    # files to write in
    fidOUT <- file(
      description = ImagePath, open = "r+b", blocking = TRUE,
      encoding = getOption("encoding"), raw = FALSE
    )
    if (!SeqRead_Image$ReadByte_Start[i] == 1) {
      nbSkip <- (SeqRead_Image$ReadByte_Start[i] - 1) * Image_Format$Bytes
      seek(fidOUT, where = nbSkip, origin = "start", rw = "write")
    }
    if (is.na(dim(ImgWrite)[3])){
      ImgChunk <- array(ImgWrite[Line_Start:(nbLines+Line_Start-1),], c(nbLines, HDR$samples, HDR$bands))
    } else {
      ImgChunk <- array(ImgWrite[Line_Start:(nbLines+Line_Start-1),,], c(nbLines, HDR$samples, HDR$bands))
    }
    ImgChunk <- aperm(ImgChunk, c(2, 3, 1))
    # writeBin(as.numeric(c(PCA_Chunk)), fidOUT, size = PCA_Format$Bytes,endian = .Platform$endian)
    if (!Image_Format$Bytes==1){
      writeBin(c(ImgChunk), fidOUT, size = Image_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    } else {
      writeBin(c(as.integer(ImgChunk)), fidOUT, size = Image_Format$Bytes, endian = .Platform$endian, useBytes = FALSE)
    }
    close(fidOUT)
  }
  rm(ImgWrite)
  rm(ImgChunk)
  gc()
  return("")
}

#' write an image resulting from "window processing" at native spatial resolution
#' (assuming square windows & origin at top left corner)
#'
#' @param Image numeric. Image corresponding to a 2D matrix or 3D array
#' @param ImagePath character. Path where the image should be written
#' @param HDR list. Image header
#' @param window_size numeric. window size used to generate Image
#'
#' @return None
#' @export

Write_Image_NativeRes <- function(Image,ImagePath,HDR,window_size){
  # update HDR: dimensions & spatial resolution
  HDR_Full <- HDR
  HDR_Full$samples <- HDR$samples * window_size
  HDR_Full$lines <- HDR$lines * window_size
  HDR_Full <- revert_resolution_HDR(HDR_Full, window_size)
  # define name for native resolution image and corresponding HDR
  ImagePath_FullRes <- paste(ImagePath, "_Fullres", sep = "")
  headerFpath <- paste(ImagePath_FullRes, ".hdr", sep = "")
  # write header
  write_ENVI_header(HDR_Full, headerFpath)
  Image_Format <- ENVI_type2bytes(HDR_Full)
  # create image the same size as native image: convert matrix into array
  if (length(dim(Image))==2){
    Image = array(Image,c(HDR$lines,HDR$samples,1))
  }
  Image_FullRes <- array(NA,c(HDR_Full$lines,HDR_Full$samples,HDR_Full$bands))
  for (band in  1:HDR_Full$bands){
    for (i in 1:HDR$lines) {
      for (j in 1:HDR$samples) {
        Image_FullRes[((i-1)*window_size+1):(i*window_size),((j-1)*window_size+1):(j*window_size),band] <- Image[i,j,band]
      }
    }
  }
  # write image and make sure size does not matter ...
  Write_Big_Image(ImgWrite = Image_FullRes,ImagePath = ImagePath_FullRes,
                  HDR = HDR_Full,Image_Format = Image_Format)
  # zip resulting file
  ZipFile(ImagePath_FullRes)
  return("")
}


#' Writes a matrix or an array into a ENVI BIL raster
#'
#' @param Image numeric. matrix or array of image to be written
#' @param HDR hdr template
#' @param ImagePath path of image file to be written
#' @param window_size spatial units use dto compute diversiy (in pixels)
#' @param FullRes should full resolution image be written (original pixel size)
#' @param LowRes should low resolution image be written (one value per spatial unit)
#' @param SmoothImage boolean. set TRUE if you want smooting filter applied to resulting diversity rasters
#
#' @return None

write_raster <- function(Image, HDR, ImagePath, window_size, FullRes = TRUE, LowRes = FALSE,SmoothImage = FALSE) {

  # check image format
  HDR$lines <- dim(Image)[1]
  HDR$samples <- dim(Image)[2]

  Image_Format <- ENVI_type2bytes(HDR)
  # Write image with resolution corresponding to window_size
  if (LowRes == TRUE) {
    # write header
    headerFpath <- paste(ImagePath, ".hdr", sep = "")
    write_ENVI_header(HDR, headerFpath)
    # write image
    Write_Big_Image(Image,ImagePath,HDR,Image_Format)
  }

  # Write image with Full native resolution
  if (FullRes == TRUE) {
    Write_Image_NativeRes(Image,ImagePath,HDR,window_size)
  }

  # write smoothed map
  if (SmoothImage == TRUE){
    SizeFilt <- 1
    if (min(c(dim(Image)[1], dim(Image)[2])) > (2 * SizeFilt + 1)) {
      Image_Smooth <- mean_filter(Image, SizeFilt,NA_remove = TRUE)
      ImagePath.Smooth <- paste(ImagePath, "_MeanFilter", sep = "")
      if (LowRes == TRUE) {
        # write header
        headerFpath <- paste(ImagePath.Smooth, ".hdr", sep = "")
        write_ENVI_header(HDR, headerFpath)
        # write image
        Write_Big_Image(Image_Smooth,ImagePath.Smooth,HDR,Image_Format)
        # close(fidOUT)
      }
      if (FullRes == TRUE) {
        # Write image with Full native resolution
        Write_Image_NativeRes(Image_Smooth,ImagePath.Smooth,HDR,window_size)
      }
    }
  }
  return("")
}


#' convert image coordinates from X-Y to index
#'
#' @param HDR_Raster list. Header file
#' @param Pixels list. row/col coordinates corresponding to the raster
#'
#' @return Image_Index
#' @export

sub2ind <- function(HDR_Raster, Pixels) {
  Image_Index <- (Pixels$col - 1) * HDR_Raster$lines + Pixels$row
  return(Image_Index)
}

#' defines the number of pieces resulting from image split
#'
#' @param HDR list. Header file
#' @param LimitSizeGb numeric. maximum size of individual pieces of an image (in Gb)
#'
#' @return nbPieces number of pieces
#' @export

split_image <- function(HDR, LimitSizeGb = FALSE) {
  Image_Format <- ENVI_type2bytes(HDR)
  lenTot <- as.double(HDR$samples) * as.double(HDR$lines) * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image_Format$Bytes) / (1024^3)
  # maximum image size read at once. If image larger, then reads in multiple pieces
  if (LimitSizeGb == FALSE) {
    LimitSizeGb <- 0.25
  }
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
  return(nbPieces)
}

#' revert resolution in a HDR file
#'
#' @param HDR information read from a header file
#' @param window_size numeric. multiplying factor for initial resolution
#'
#' @return updated HDR information
#' @export
revert_resolution_HDR <- function(HDR, window_size) {
  MapInfo <- strsplit(HDR$`map info`, split = ",")
  MapInfo[[1]][6] <- as.numeric(MapInfo[[1]][6]) / window_size
  MapInfo[[1]][7] <- as.numeric(MapInfo[[1]][7]) / window_size
  HDR$`map info` <- paste(MapInfo[[1]], collapse = ",")
  return(HDR)
}

#' Zips an image file
#'
#' @param ImagePath character. path for the image
#' @return None
#' @export
ZipFile <- function(ImagePath) {


  ImagePath <- normalizePath(ImagePath)

  zip::zipr(zipfile = paste0(ImagePath, ".zip"), files = ImagePath)
  file.remove(ImagePath)

  return(invisible())
}
