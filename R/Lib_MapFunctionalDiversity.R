# ==============================================================================
# biodivMapR
# Lib_MapFunctionalDiversity.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This Library produces maps corresponding to functional diversity indicators
# (Richness, Evenness & Divergence), following the definition of
# Villeger, Mason & Mouillot (2008), NEW MULTIDIMENSIONAL FUNCTIONAL DIVERSITY
# INDICES FOR A MULTIFACETED FRAMEWORK IN FUNCTIONAL ECOLOGY, Ecology
# (https://doi.org/10.1890/07-1206.1)
# based on a set of descriptors produced from a raster with biodivMapR (usually selected components from PCA)
# ==============================================================================

#' maps functional diversity indicators based on prior selection of PCs
#'
#' @param Original_Image_File character. Path and name of the original input image for biodivMapR.
#' @param Functional_File character. Path and name of the image processed to be used to compute functional diversity
#' --> can be PCA file with. if FALSE: using Original_Image_File
#' @param Selected_Features numeric. Contains features to be used from Input_Image_File. using all if FALSE
#' @param Output_Dir character. Output directory.
#' @param window_size numeric. Size of spatial units (in pixels) to compute diversity.
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...).
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param FullRes boolean. Full resolution.
#' @param LowRes boolean. Low resolution.
#' @param MapSTD boolean. map of standard deviation of the alpha diversity map (over repetitions)
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param SmoothImage boolean. set TRUE if you want smooting filter applied to resulting diversity rasters
#' @param FDmetric character. Functional diversity metric
#'
#' @return None
#' @importFrom future plan multisession sequential
#' @export
#'
map_functional_div <- function(Original_Image_File,Functional_File = FALSE,
                               Selected_Features = FALSE,
                               Output_Dir, window_size,
                               TypePCA = "SPCA",
                               MinSun = 0.25,
                               FullRes = FALSE, LowRes = TRUE, MapSTD = TRUE,
                               nbCPU = 1, MaxRAM = 0.25,
                               SmoothImage = TRUE,
                               FDmetric = c('FRic', 'FEve', 'FDiv')) {

  if (Functional_File==FALSE) Functional_File <- Original_Image_File
  # check if selected features match with image dimensions
  HDRname <- get_HDR_name(Functional_File)
  HDR <- read_ENVI_header(HDRname)
  if (length(Selected_Features) ==1 & Selected_Features[1] ==FALSE){
    Selected_Features = seq(1,HDR$bands)
  } else {
    if (max(Selected_Features)>HDR$bands){
      message("*********************************************************")
      message("  WARNING: Selected_Features includes more features than ")
      message("                 available in input file                 ")
      print(Functional_File)
      message("                     process aborted                     ")
      message("*********************************************************")
      stop()
    }
  }
  # define output directory
  Output_Dir_Funct <- define_output_subdir(Output_Dir, Original_Image_File, TypePCA, "FUNCTIONAL")
  # 1- COMPUTE FUNCTIONAL DIVERSITY: RICHNESS, EVENNESS, DIVERGENCE
  print("Compute functional metrics")
  FunctionalMetrics <- compute_Functional_metrics(Functional_File = Functional_File,
                                                  Functional_Map_Path = Functional_Map_Path,
                                                  Selected_Features = Selected_Features,
                                                  FDmetric = FDmetric,
                                                  window_size = window_size,
                                                  MinSun = MinSun, nbCPU = nbCPU, MaxRAM = MaxRAM)

  ## prepare header for functional diversity map
  HDRname <- get_HDR_name(Functional_File)
  HDR_Funct <- read_ENVI_header(HDRname)
  HDR_Funct$bands <- 1
  HDR_Funct$`data type` <- 4
  # define image size
  HDR_Funct$lines <- dim(FunctionalMetrics[[1]])[1]
  HDR_Funct$samples <- dim(FunctionalMetrics[[1]])[2]
  # change resolution
  HDR_Funct <- change_resolution_HDR(HDR_Funct, window_size)
  # 2- SAVE FUNCTIONAL DIVERSITY MAPS
  print("Write functional diversity maps")
  for (FD in FDmetric){
    Functional_Map_Path <- file.path(Output_Dir_Funct, FD)
    HDR_Funct$`band names` <- FD
    write_raster(Image = FunctionalMetrics[[FD]],
                 HDR = HDR_Funct,
                 ImagePath = Functional_Map_Path,
                 window_size = window_size,
                 FullRes = FullRes,
                 LowRes = LowRes,
                 SmoothImage = SmoothImage)
  }
  return(invisible())
}

#' Map functional diversity metrics based on spectral species
#'
#' @param Functional_File character. Path and name of the image processed to be used to compute functional diversity
#' @param Functional_Map_Path character. Path and name of the resulting functional diversity raster
#' @param Selected_Features numeric. list of features from Functional_File to be included in the analysis
#' @param FDmetric character. Functional diversity metric
#' @param window_size numeric. size of spatial units (in pixels) to compute diversity
#' @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#
#' @return list of mean and SD of alpha diversity metrics
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats sd
#' @import cli
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom terra rast values
#' @export

compute_Functional_metrics <- function(Functional_File, Functional_Map_Path, Selected_Features,
                                       FDmetric = c('FRic', 'FEve', 'FDiv'),
                                       window_size, MinSun, nbCPU = 1, MaxRAM = 0.25) {

  ## read Functional_File and write Functional_Map_Path
  HDRname <- get_HDR_name(Functional_File)
  HDR <- read_ENVI_header(HDRname)
  nbPieces_Min <- split_image(HDR, MaxRAM)
  if (nbPieces_Min < nbCPU) nbPieces_Min <- nbCPU
  SeqRead.Funct <- where_to_read_kernel(HDR, nbPieces_Min, window_size)

  # for each piece of image
  ReadWrite <- list()
  for (i in 1:nbPieces_Min) {
    ReadWrite[[i]] <- list()
    ReadWrite[[i]]$Byte_Start <- SeqRead.Funct$ReadByte_Start[i]
    ReadWrite[[i]]$nbLines <- SeqRead.Funct$Lines_Per_Chunk[i]
    ReadWrite[[i]]$lenBin <- SeqRead.Funct$ReadByte_End[i] - SeqRead.Funct$ReadByte_Start[i] + 1
  }
  ImgFormat <- "3D"
  FunctIN_Format <- ENVI_type2bytes(HDR)

  # get interquantile range for standardization
  MinFunct <- MaxFunct <- c()
  for (i in Selected_Features){
    RasterVal <- terra::rast(Functional_File, lyrs  = i)
    RangeTraits <- IQR_outliers(terra::values(RasterVal))
    MinFunct <- c(MinFunct, RangeTraits[1])
    MaxFunct <- c(MaxFunct, RangeTraits[2])
  }
  MinMaxRaster <- data.frame('MinRaster'=MinFunct,'MaxRaster'=MaxFunct)

  message(paste('compute Functional diversity in', nbPieces_Min, 'chunks'))
  FUNCT_DIV <- lapply(ReadWrite, FUN = Get_FunctionalMetrics_From_Traits,
                      Functional_File = Functional_File,
                      Selected_Features = Selected_Features,
                      MinMaxRaster = MinMaxRaster,
                      HDR = HDR, FunctIN_Format = FunctIN_Format,
                      ImgFormat = ImgFormat,
                      window_size = window_size,
                      MinSun = MinSun,
                      FDmetric = FDmetric,
                      nbCPU = nbCPU)

  # create ful map from chunks
  FRic_Chunk <- FEve_Chunk <- FDiv_Chunk <- list()
  for (i in 1:length(FUNCT_DIV)) {
    FRic_Chunk[[i]] <- FUNCT_DIV[[i]]$FRic
    FEve_Chunk[[i]] <- FUNCT_DIV[[i]]$FEve
    FDiv_Chunk[[i]] <- FUNCT_DIV[[i]]$FDiv
  }
  FRicMap <- do.call(rbind, FRic_Chunk)
  FEveMap <- do.call(rbind, FEve_Chunk)
  FDivMap <- do.call(rbind, FDiv_Chunk)
  FunctMap <- list()
  FunctMap$FRic <- FRicMap
  FunctMap$FEve <- FEveMap
  FunctMap$FDiv <- FDivMap
  return(FunctMap)
  # return(list("FunctMap" = FunctMap, "HDR" = HDR_Funct))
}

#' Prepare for the computation of the functional diversity metrics
#'
#' @param ReadWrite numeric. bytes coordinates for each read and write
#' @param Functional_File character. path for the raster file to get functional traits from
#' @param Selected_Features numeric. features to be used from Functional_File
#' @param MinMaxRaster numeric. min and max values to be used for data standardization
#' @param HDR list. header file
#' @param FunctIN_Format list. image format (bytes, data type, etc)
#' @param ImgFormat character. define if image is 2D or 3D
#' @param window_size numeric. window size to compute metrics from
#' @param MinSun numeric. minimum sun exposition to keep window as valid data (ratio of sunlit pixels / total pixels in the window)
#' @param FDmetric character. enumerate diversity metrics to be computed
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#'
#' @return FDmetrics
#' @export

Get_FunctionalMetrics_From_Traits <- function(ReadWrite, Functional_File, Selected_Features,
                                              MinMaxRaster, HDR,
                                              FunctIN_Format,
                                              ImgFormat, window_size,
                                              MinSun,
                                              FDmetric = c('FRic', 'FEve', 'FDiv'),
                                              nbCPU = 1) {

  Image_Chunk <- read_BIL_image_subset(Functional_File, HDR,
                                       ReadWrite$Byte_Start, ReadWrite$lenBin,
                                       ReadWrite$nbLines, FunctIN_Format, ImgFormat)
  Image_Chunk <- Image_Chunk[,,Selected_Features]
  if (length(Selected_Features)==1) Image_Chunk <- array(Image_Chunk,
                                                         dim = c(nrow(Image_Chunk),
                                                                 ncol(Image_Chunk),
                                                                 1))
  # standardize data
  for (i in 1:length(Selected_Features)){
    Image_Chunk[,,i] <- (Image_Chunk[,,i]-MinMaxRaster$MinRaster[i])/(MinMaxRaster$MaxRaster[i]-MinMaxRaster$MinRaster[i])
  }
  FDmetrics <- compute_FD(Image_Chunk = Image_Chunk,
                          window_size = window_size,
                          MinSun = MinSun,
                          FDmetric = FDmetric,
                          nbCPU = nbCPU)
  rm(Image_Chunk)
  gc()
  return(FDmetrics)
}

#' compute functional diversity metrics for an array, given a specific window size
#'
#' @param Image_Chunk numeric. 3D image chunk of spectral species
#' @param window_size numeric. size of spatial units (in pixels) to compute diversity
#' @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#' @param FDmetric character. Functional diversity metric
#' @param nbCPU numeric. nb of CPU to process data
#'
#' @return list of functional diversity metrics corresponding to image chunk
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom stats na.omit
#' @export

compute_FD <- function(Image_Chunk, window_size, MinSun,
                       FDmetric = c('FRic', 'FEve', 'FDiv'), nbCPU) {
  nbi <- floor(dim(Image_Chunk)[1] / window_size)
  nbj <- floor(dim(Image_Chunk)[2] / window_size)
  nbTraits <- dim(Image_Chunk)[3]
  FRicmap <- FEvemap <- FDivmap <- matrix(NA, nrow = nbi, ncol = nbj)

  listWindows <- listCoords <- list()
  ij <- 0
  # for each kernel in the line
  for (ii in 1:nbi) {
    # for each kernel in the column
    for (jj in 1:nbj) {
      li <- ((ii - 1) * window_size) + 1
      ui <- ii * window_size
      lj <- ((jj - 1) * window_size) + 1
      uj <- jj * window_size

      wintmp <- data.frame(matrix(Image_Chunk[li:ui, lj:uj, ], ncol = nbTraits))
      wintmp <- stats::na.omit(wintmp)
      nbPix_Sunlit <- dim(wintmp)[1]
      PCsun <- nbPix_Sunlit / window_size**2
      if (PCsun > MinSun & nbPix_Sunlit > 3) {
        # convert into 2D matrix shape
        ij <- ij +1
        # keep non zero values
        row.names(wintmp) <- paste('pix#',seq(1,nrow(wintmp)),sep = '')
        listWindows[[ij]] <- wintmp
        listCoords[[ij]] <- list('Row' = ii, 'Col' = jj)
      }
    }
  }
  plan(multisession, workers = nbCPU)
  handlers(global = TRUE)
  handlers("cli")
  with_progress({
    p <- progressr::progressor(steps = ij)
    alphaSSD <- future.apply::future_lapply(X = listWindows,
                                            FUN = getFD, p = p)
  })
  plan(sequential)
  for (i in 1:length(listCoords)){
    FRicmap[listCoords[[i]]$Row, listCoords[[i]]$Col] <- alphaSSD[[i]]$FRic
    FEvemap[listCoords[[i]]$Row, listCoords[[i]]$Col] <- alphaSSD[[i]]$FEve
    FDivmap[listCoords[[i]]$Row, listCoords[[i]]$Col] <- alphaSSD[[i]]$FDiv
  }
  my_list <- list("FRic" = FRicmap, "FDiv" = FDivmap, "FEve" = FEvemap)
  return(my_list)
}


#' get functional diversity metrics from dataframe
#' This function was inspired from FD package
#' @param spectraits numeric. dataframe containing species in rows and trait values in columns
#' @param FDmetric character. Functional diversity metric
#' @param p list. progressor object for progress bar
#
#' @return FDmetrics
#' @importFrom geometry convhulln
#' @importFrom ape mst
#' @importFrom stats dist
#'
#' @export

getFD <- function(spectraits,
                  FDmetric = c('FRic', 'FEve', 'FDiv'),
                  p = NULL){
  nbTraits <- ncol(spectraits)
  nbSpecies <- nrow(spectraits)
  FRic <- FDiv <- FEve <- NA
  if (nbTraits>1){
    if (!is.na(match('FRic', FDmetric)) | !is.na(match('FDiv', FDmetric))){
      # inspired from FD package
      # convex hull using geometry
      FunctHull <- geometry::convhulln(spectraits, "Fx TO 'vert.txt'", output.options = 'FA')
      # 1- Functional Richness
      if (!is.na(match('FRic', FDmetric))){
        FRic <- FunctHull$vol
      }
      # 2- Functional Divergence mean distance from centroid
      if (!is.na(match('FDiv', FDmetric))){
        vert1 <- scan("vert.txt", quiet = T)
        vert2 <- vert1 + 1
        vertices <- vert2[-1]
        trvertices <- spectraits[vertices, ]
        # coordinates of the center of gravity of the vertices (Gv)
        baryv <- apply(trvertices, 2, mean)
        # euclidian dstances to Gv (dB)
        distbaryv <- rep(0, nbSpecies)
        for (j in 1:nbSpecies) distbaryv[j] <- (sum((spectraits[j, ] - baryv)^2) ) ^0.5
        # mean of dB values
        meandB <- mean(distbaryv)
        # deviations to mean of db
        devdB <- distbaryv - meandB
        # computation of FDiv
        FDiv <- (sum(devdB/nbSpecies) + meandB) / (sum(abs(devdB/nbSpecies)) + meandB)
      }
    }
    if (!is.na(match('FEve', FDmetric))){
      # computation of minimum spanning tree and conversion of the 'mst' matrix into 'dist' class
      tr.dist <- stats::dist(spectraits)
      linkmst <- ape::mst(tr.dist)
      mstvect <- as.dist(linkmst)
      # computation of EW for the (nbSpecies - 1) segments to link the nbSpecies points
      EW <- rep(0, nbSpecies - 1)
      flag <- 1
      for (m in 1 : ((nbSpecies - 1) * nbSpecies / 2)) {
        if (mstvect[m] != 0) {
          EW[flag] <- tr.dist[m]
          flag <- flag + 1
        }
      }
      # computation of the PEW and comparison with 1 / nbSpecies - 1, finally computation of FEve
      minPEW <- rep(0, nbSpecies - 1)
      OdSmO <- 1 / (nbSpecies - 1)
      for (l in 1 : (nbSpecies - 1)) {
        minPEW[l] <- min((EW[l] / sum(EW)), OdSmO)
      }
      FEve <- ((sum(minPEW)) - OdSmO) / (1 - OdSmO)
    }
  }
  if (!is.null(p)){p()}
  return(list('FRic' = FRic, 'FDiv' = FDiv, 'FEve' = FEve))
}
