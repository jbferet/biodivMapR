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
#'
#' @return None
#' @importFrom future plan multiprocess sequential
#' @export
#'
map_functional_div <- function(Original_Image_File,Functional_File = FALSE,
                               Selected_Features = FALSE,
                               Output_Dir, window_size,
                               TypePCA = "SPCA",
                               MinSun = 0.25,
                               FullRes = TRUE,LowRes = FALSE, MapSTD = FALSE,
                               nbCPU = FALSE, MaxRAM = FALSE,SmoothImage = TRUE) {

  if (Functional_File==FALSE){
    Functional_File <- Original_Image_File
  }
  # check if selected features match with image dimensions
  HDRname <- get_HDR_name(Functional_File)
  HDR <- read_ENVI_header(HDRname)
  if (length(Selected_Features) ==1 & Selected_Features[1] ==FALSE){
    Selected_Features = seq(1,HDR$bands)
  } else {
    if (max(Selected_Features)>HDR$bands){
      message("")
      message("*********************************************************")
      message("  WARNING: Selected_Features includes more features than ")
      message("                 available in input file                 ")
      print(Functional_File)
      message("process aborted")
      message("*********************************************************")
      message("")
      stop()
    }
  }
  # define output directory
  Output_Dir_Funct <- define_output_subdir(Output_Dir, Original_Image_File, TypePCA, "FUNCTIONAL")
  Functional_Map_Path <- file.path(Output_Dir_Funct, "FunctionalDiversity_Map")
  # 1- COMPUTE FUNCTIONAL DIVERSITY: RICHNESS, EVENNESS, DIVERGENCE
  print("Compute functional metrics")
  FunctionalMetrics <- compute_Functional_metrics(Functional_File = Functional_File, Functional_Map_Path = Functional_Map_Path,
                                                  Selected_Features = Selected_Features, window_size = window_size,
                                                  MinSun = MinSun, nbCPU = nbCPU, MaxRAM = MaxRAM)

  # 2- SAVE FUNCTIONAL DIVERSITY MAPS
  print("Write functional diversity maps")
  write_raster(FunctionalMetrics$FunctMap, FunctionalMetrics$HDR, Functional_Map_Path,
               window_size, FullRes = FullRes, LowRes = LowRes,SmoothImage = SmoothImage)
  return(invisible())
}

#' Map functional diversity metrics based on spectral species
#'
#' @param Functional_File character. Path and name of the image processed to be used to compute functional diversity
#' @param Functional_Map_Path character. Path and name of the resulting functional diversity raster
#' @param Selected_Features numeric. list of features from Functional_File to be included in the analysis
#' @param window_size numeric. size of spatial units (in pixels) to compute diversity
#' @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#
#' @return list of mean and SD of alpha diversity metrics
#' @importFrom future plan multiprocess sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats sd
#' @importFrom raster brick values nbands
compute_Functional_metrics <- function(Functional_File, Functional_Map_Path, Selected_Features,
                                       window_size, MinSun, nbCPU = FALSE, MaxRAM = FALSE) {

  ## read Functional_File and write Functional_Map_Path
  HDRname <- get_HDR_name(Functional_File)
  HDR <- read_ENVI_header(HDRname)
  if (MaxRAM == FALSE) {
    MaxRAM <- 0.25
  }
  nbPieces_Min <- split_image(HDR, MaxRAM)
  if (nbCPU == FALSE) {
    nbCPU <- 1
  }
  if (nbPieces_Min < nbCPU) {
    nbPieces_Min <- nbCPU
  }
  SeqRead.Funct <- where_to_read_kernel(HDR, nbPieces_Min, window_size)

  ## prepare functional diversity map
  HDR_Funct <- HDR
  # define number of bands for functional diversity
  HDR_Funct$bands <- 3
  HDR_Funct$`data type` <- 4
  # define image size
  HDR_Funct$samples <- floor(HDR_Funct$samples / window_size)
  HDR_Funct$lines <- floor(HDR_Funct$lines / window_size)
  # change resolution
  HDR_Funct <- change_resolution_HDR(HDR_Funct, window_size)
  HDR_Funct$`band names` <- c('RICHNESS','EVENNESS','DIVERGENCE')
  SeqWrite.FUNCT <- where_to_write_kernel(HDR, HDR_Funct, nbPieces_Min, window_size)

  # for each piece of image
  ReadWrite <- list()
  for (i in 1:nbPieces_Min) {
    ReadWrite[[i]] <- list()
    ReadWrite[[i]]$RW_FUNCTmap <- ReadWrite[[i]]$RW_FUNCT <- list()
    ReadWrite[[i]]$RW_FUNCT$Byte_Start <- SeqRead.Funct$ReadByte_Start[i]
    ReadWrite[[i]]$RW_FUNCT$nbLines <- SeqRead.Funct$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW_FUNCT$lenBin <- SeqRead.Funct$ReadByte_End[i] - SeqRead.Funct$ReadByte_Start[i] + 1

    ReadWrite[[i]]$RW_FUNCTmap$Byte_Start <- SeqWrite.FUNCT$ReadByte_Start[i]
    ReadWrite[[i]]$RW_FUNCTmap$nbLines <- SeqWrite.FUNCT$Lines_Per_Chunk[i]
    ReadWrite[[i]]$RW_FUNCTmap$lenBin <- SeqWrite.FUNCT$ReadByte_End[i] - SeqWrite.FUNCT$ReadByte_Start[i] + 1
  }
  ImgFormat <- "3D"
  FunctOUT_Format <- ENVI_type2bytes(HDR_Funct)
  FunctIN_Format <- ENVI_type2bytes(HDR)

  # # get minimum and maximum value for each feature
  # FunctRaster <- brick(Functional_File)
  # RasterVal <- setMinMax(FunctRaster)
  # MinFunct <- MaxFunct <- c()
  # for (i in 1:nbands(FunctRaster)){
  #   MinFunct[i] <- minValue(RasterVal[[i]])
  #   MaxFunct[i] <- maxValue(RasterVal[[i]])
  # }

  # get interquantile range for standardization
  FunctRaster <- brick(Functional_File)
  MinFunct <- MaxFunct <- c()
  for (i in 1:nbands(FunctRaster)){
    RasterVal <- FunctRaster[[i]]
    RangeTraits <- IQR_outliers(raster::values(RasterVal))
    MinFunct[i] <- RangeTraits[1]
    MaxFunct[i] <- RangeTraits[2]
  }

  MinFunct <- MinFunct[Selected_Features]
  MaxFunct <- MaxFunct[Selected_Features]
  MinMaxRaster <- data.frame('MinRaster'=MinFunct,'MaxRaster'=MaxFunct)

  # multiprocess of spectral species distribution and alpha diversity metrics
  if (nbCPU>1){
    plan(multiprocess, workers = nbCPU) ## Parallelize using four cores
    Schedule_Per_Thread <- ceiling(nbPieces_Min / nbCPU)
    FUNCT_DIV <- future_lapply(ReadWrite,
                               FUN = Get_FunctionalMetrics_From_Traits, Functional_File = Functional_File, Selected_Features = Selected_Features,
                               MinMaxRaster = MinMaxRaster, HDR = HDR, HDR_Funct = HDR_Funct,
                               FunctIN_Format = FunctIN_Format, FunctOUT_Format = FunctOUT_Format,
                               ImgFormat = ImgFormat, window_size = window_size, MinSun = MinSun,
                               Functional_Map_Path = Functional_Map_Path,  future.scheduling = Schedule_Per_Thread
    )
    plan(sequential)
  } else {
    FUNCT_DIV <- lapply(ReadWrite, FUN = Get_FunctionalMetrics_From_Traits, Functional_File = Functional_File,
                        Selected_Features = Selected_Features, MinMaxRaster = MinMaxRaster,
                        HDR = HDR, HDR_Funct = HDR_Funct, FunctIN_Format = FunctIN_Format,
                        FunctOUT_Format = FunctOUT_Format, ImgFormat = ImgFormat, window_size = window_size,
                        MinSun = MinSun, Functional_Map_Path = Functional_Map_Path)
  }
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
  HDR_Funct$lines <- dim(FRicMap)[1]
  HDR_Funct$samples <- dim(FRicMap)[2]
  FunctMap <- array(NA, c(HDR_Funct$lines, HDR_Funct$samples, 3))
  FunctMap[,,1] <- FRicMap
  FunctMap[,,2] <- FEveMap
  FunctMap[,,3] <- FDivMap
  my_list <- list("FunctMap" = FunctMap, "HDR" = HDR_Funct)
  return(my_list)
}

# Prepare for the computation of the functional diversity metrics
#
# @param ReadWrite
# @param Functional_File
# @param Selected_Features
# @param MinMaxRaster
# @param HDR
# @param HDR_Funct
# @param FunctIN_Format
# @param FunctOUT_Format
# @param ImgFormat
# @param window_size
# @param MinSun
# @param Functional_Map_Path
#
# @importFrom raster brick nbands minValue maxValue
# @return
Get_FunctionalMetrics_From_Traits <- function(ReadWrite, Functional_File, Selected_Features,
                                              MinMaxRaster, HDR, HDR_Funct,
                                              FunctIN_Format, FunctOUT_Format, ImgFormat, window_size,
                                              MinSun, Functional_Map_Path) {

  Image_Chunk <- read_BIL_image_subset(Functional_File, HDR,
                                   ReadWrite$RW_FUNCT$Byte_Start, ReadWrite$RW_FUNCT$lenBin,
                                   ReadWrite$RW_FUNCT$nbLines, FunctIN_Format, ImgFormat)
  Image_Chunk <- Image_Chunk[,,Selected_Features]
  # standardize data
  for (i in 1:length(Selected_Features)){
    Image_Chunk[,,i] <- (Image_Chunk[,,i]-MinMaxRaster$MinRaster[i])/(MinMaxRaster$MaxRaster[i]-MinMaxRaster$MinRaster[i])
  }
  FUNCTmetrics <- compute_FUNCT(Image_Chunk, window_size, MinSun)
  rm(Image_Chunk)
  gc()
  return(FUNCTmetrics)
}

#' compute functional diversity metrics for an array, given a specific window size
#'
#' @param Image_Chunk numeric. 3D image chunk of spectral species
#' @param window_size numeric. size of spatial units (in pixels) to compute diversity
#' @param MinSun numeric. minimum proportion of sunlit pixels required to consider plot
#'
#' @return list of functional diversity metrics corresponding to image chunk
#' @importFrom geometry convhulln
# @importFrom emstreeR ComputeMST
#' @export

compute_FUNCT <- function(Image_Chunk, window_size, MinSun) {
  nbi <- floor(dim(Image_Chunk)[1] / window_size)
  nbj <- floor(dim(Image_Chunk)[2] / window_size)
  nbTraits <- dim(Image_Chunk)[3]
  FRicmap <- FEvemap <- FDivmap <- matrix(NA, nrow = nbi, ncol = nbj)

  # for each kernel in the line
  for (ii in 1:nbi) {
    # for each kernel in the column
    for (jj in 1:nbj) {
      li <- ((ii - 1) * window_size) + 1
      ui <- ii * window_size
      lj <- ((jj - 1) * window_size) + 1
      uj <- jj * window_size
      # convert into 2D matrix shape
      ij <- matrix(Image_Chunk[li:ui, lj:uj, ], ncol = nbTraits)
      # keep non zero values
      ij <- matrix(ij[which(!is.na(ij[,1])),], ncol = nbTraits)
      nbPix_Sunlit <- dim(ij)[1]
      PCsun <- nbPix_Sunlit / window_size**2
      if (PCsun > MinSun) {
        # compute functional metrics
        # 1- Functional Richness
        # convex hull using geometry
        FRicmap[ii,jj] <- 100*geometry::convhulln(ij, output.options = 'FA')$vol
        # 2- Functional Divergence
        # mean distance from centroid
        Centroid <- colMeans(ij)
        FDivmap[ii,jj] <- 100*sum(sqrt(rowSums((t(t(ij) - Centroid)^2))))/nbPix_Sunlit
        # FDivmap[ii,jj] <- 100*sum(sqrt(rowSums((t(t(ij) )^2))))/nbPix_Sunlit
        # 3- Functional Evenness
        # euclidean minimum spanning tree
        # 20220112: wait for update of emstreeR and integration of mlpack
        # FEvemap[ii,jj] <- 100*sum(emstreeR::ComputeMST(ij,verbose = FALSE)$distance)/nbPix_Sunlit
        FEvemap[ii,jj] <- NA*FDivmap[ii,jj]
      } else {
        FRicmap[ii,jj] <- NA
        FDivmap[ii,jj] <- NA
        FEvemap[ii,jj] <- NA
      }
    }
  }
  my_list <- list("FRic" = FRicmap, "FDiv" = FDivmap, "FEve" = FEvemap)
  return(my_list)
}
