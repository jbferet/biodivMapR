# ==============================================================================
# biodivMapR
# Lib_MapPartition.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2020/05 Jean-Baptiste FERET
# ==============================================================================
# This Library produces maps corresponding to the Partitioning of plant spectral
# diversity into alpha and beta components, following the method proposed by
# Laliberte, Schweiger & Legendre (2020), Partitioning plant spectral diversity
# into alpha and beta components, Ecology letters (https://doi.org/10.1111/ele.13429)
# ==============================================================================

#' Partitioning of plant spectral diversity
#'
#' @param Original_Image_File character. Path and name of the original input image for biodivMapR.
#' @param Partition_File character. Path and name of the input image to perform diversity partitioning on
#' @param Selected_Features numeric. Contains features to be used from Input_Image_File. using all if FALSE
#' @param Output_Dir character. Output directory.
#' @param window_size numeric. Size of spatial units (in pixels) to compute diversity.
#' @param TypePCA character. Type of PCA (PCA, SPCA, NLPCA...).
#' @param MinSun numeric. Minimum proportion of sunlit pixels required to consider plot.
#' @param FullRes boolean. Full resolution.
#' @param LowRes boolean. Low resolution.
#' @param nbCPU numeric. Number of CPUs to use in parallel.
#' @param MaxRAM numeric. MaxRAM maximum size of chunk in GB to limit RAM allocation when reading image file.
#' @param SmoothImage boolean. set TRUE if you want smooting filter applied to resulting diversity rasters
#'
#' @return None
#' @importFrom raster brick
#' @export
#'
map_partition_div <- function(Original_Image_File,Partition_File = FALSE,
                              Selected_Features = FALSE,Output_Dir, window_size,
                              TypePCA = "SPCA", MinSun = 0.25,
                              FullRes = TRUE,LowRes = FALSE, nbCPU = FALSE, MaxRAM = FALSE,
                              SmoothImage = TRUE) {

  if (Partition_File==FALSE){
    Partition_File <- Original_Image_File
  }
  # check if selected features match with image dimensions
  HDRname <- get_HDR_name(Partition_File)
  HDR <- read_ENVI_header(HDRname)
  if (Selected_Features==FALSE){
    Selected_Features = seq(1,HDR$bands)
  } else {
    if (max(Selected_Features)>HDR$bands){
      message("")
      message("*********************************************************")
      message("  WARNING: Selected_Features includes more features than ")
      message("                 available in input file                 ")
      print(Partition_File)
      message("process aborted")
      message("*********************************************************")
      message("")
      stop()
    }
  }
  # define output directory
  Output_Dir_Partition <- define_output_subdir(Output_Dir, Original_Image_File, TypePCA, "PARTITION_SPECTRAL")
  print("Partitioning spectral diversity following Lalibert?, Schweiger & Legendre (2020)")
  brick_SelPC <- brick(Partition_File)
  if (!Selected_Features == FALSE){
    brick_SelPC <- brick(brick_SelPC[[Selected_Features]])
  }
  len <- window_size*window_size
  Partition_Spectral <- specdiv(brick_SelPC, fact = window_size, prop = MinSun)

  print("Writing alpha spectral diversity")
  ## prepare HDR
  HDR_Partition <- HDR
  HDR_Partition$bands <- 1
  HDR_Partition$`data type` <- 4
  HDR_Partition$samples <- dim(Partition_Spectral$rasters$alpha_sdiv)[2]
  HDR_Partition$lines <- dim(Partition_Spectral$rasters$alpha_sdiv)[1]
  HDR_Partition <- change_resolution_HDR(HDR_Partition, window_size)
  HDR_Partition$`band names` <- c('Alpha_spectral_diversity')
  # Write image
  RasterPart <- t(matrix(raster::values(Partition_Spectral$rasters$alpha_sdiv), nrow = HDR_Partition$samples, ncol = HDR_Partition$lines))
  Alpha_Map_Path <- paste(Output_Dir_Partition, "Alpha_SpectralDiv", sep = "")
  write_raster(RasterPart, HDR_Partition, Alpha_Map_Path, window_size, FullRes = FullRes, LowRes = LowRes,SmoothImage = SmoothImage)

  print("Writing local contribution to spectral gamma-diversity")
  ## prepare HDR
  HDR_Partition$`band names` <- c('LCSD')
  # Write image
  RasterPart <- t(matrix(raster::values(Partition_Spectral$rasters$beta_lcsd), nrow = HDR_Partition$samples, ncol = HDR_Partition$lines))
  LCSD_Map_Path <- paste(Output_Dir_Partition, "LCSD", sep = "")
  write_raster(RasterPart, HDR_Partition, LCSD_Map_Path, window_size, FullRes = FullRes, LowRes = LowRes,SmoothImage = SmoothImage)

  print("Writing LCSS")
  ## prepare HDR
  HDR_Partition$`band names` <- c('LCSS')
  # Write image
  RasterPart <- t(matrix(raster::values(Partition_Spectral$rasters$beta_lcss), nrow = HDR_Partition$samples, ncol = HDR_Partition$lines))
  LCSS_Map_Path <- paste(Output_Dir_Partition, "LCSS", sep = "")
  write_raster(RasterPart, HDR_Partition, LCSS_Map_Path, window_size, FullRes = FullRes, LowRes = LowRes,SmoothImage = SmoothImage)

  print("Writing feature contribution to spectral alpha diversity")
  ## prepare HDR
  HDR_Partition$bands <- length(Selected_Features)
  HDR_Partition$`band names` <- paste('Contribution_Feature', 1:HDR_Partition$bands, collapse = ", ")
  # Write image
  RasterPart <- aperm(array(raster::values(Partition_Spectral$rasters$alpha_fcsd), c(HDR_Partition$samples, HDR_Partition$lines, HDR_Partition$bands)),c(2,1,3))
  FCSD_Map_Path <- paste(Output_Dir_Partition, "FCSD", sep = "")
  write_raster(RasterPart, HDR_Partition, FCSD_Map_Path, window_size, FullRes = FullRes, LowRes = LowRes,SmoothImage = SmoothImage)
  my_list <- list("Sum_Squares" = Partition_Spectral$ss,
                  "Spectral_Div" = Partition_Spectral$sdiv,
                  "FeatureContribution" = Partition_Spectral$fcsd)
  return(my_list)
}
