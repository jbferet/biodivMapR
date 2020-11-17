# ==============================================================================
# biodivMapR
# Lib_FilterData.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Florian de Boissieu <fdeboiss@gmail.com>
# Copyright 2020/06 Jean-Baptiste FERET
# ==============================================================================
# This library contains functions to filter raster based on radiometric criteria
# ==============================================================================

#' Performs radiometric filtering based on three criteria: NDVI, NIR reflectance, Blue reflectance
#'
#' @param Image_Path character. Path of the image to be processed
#' @param Mask_Path character. Path of the mask corresponding to the image
#' @param Output_Dir character. Path for output directory
#' @param TypePCA character. Type of PCA: choose either "PCA" or "SPCA"
#' @param NDVI_Thresh numeric. NDVI threshold applied to produce a mask (select pixels with NDVI>NDVI_Thresh)
#' @param Blue_Thresh numeric. Blue threshold applied to produce a mask (select pixels with Blue refl < Blue_Thresh --> filter clouds) refl expected between 0 and 10000
#' @param NIR_Thresh numeric. NIR threshold applied to produce a mask (select pixels with NIR refl > NIR_Thresh) refl expected between 0 and 10000
#' @param Blue numeric. central wavelength corresponding to the blue spectral band (in nanometers)
#' @param Red numeric. central wavelength corresponding to the red spectral band (in nanometers)
#' @param NIR numeric. central wavelength corresponding to the NIR spectral band (in nanometers)
#'
#' @return MaskPath = updated mask file
#' @export
perform_radiometric_filtering <- function(Image_Path, Mask_Path, Output_Dir, TypePCA = "SPCA",
                                          NDVI_Thresh = 0.5, Blue_Thresh = 500, NIR_Thresh = 1500,
                                          Blue = 480, Red = 700, NIR = 835) {
  # check if format of raster data is as expected
  check_data(Image_Path)
  if (!Mask_Path==FALSE){
    check_data(Mask_Path,Mask = TRUE)
  }
  # define full output directory
  Output_Dir_Full <- define_output_directory(Output_Dir, Image_Path, TypePCA)
  # Create / Update shade mask if optical data
  if (Mask_Path == FALSE | Mask_Path == "") {
    print("Create mask based on NDVI, NIR and Blue threshold")
  } else {
    print("Update mask based on NDVI, NIR and Blue threshold")
  }
  Shade_Update <- paste(Output_Dir_Full, "ShadeMask_Update", sep = "")
  Mask_Path <- create_mask_from_threshold(ImPath = Image_Path, MaskPath = Mask_Path, MaskPath_Update = Shade_Update,
                                          NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh, NIR_Thresh = NIR_Thresh,
                                          Blue = Blue, Red = Red, NIR = NIR)
  return(Mask_Path)
}

#' create a mask based on NDVI, Green reflectance and NIR reflectance
#' NDVI (min) threshold eliminates non vegetated pixels
#' Blue (max) threshold eliminates Clouds
#' NIR (min) threshold eliminates shadows
#' ! only valid if Optical data!!
#'
#' @param ImPath character. Full path of a raster file
#' @param MaskPath character. Full path of the mask to be used with the raster file
#' @param MaskPath_Update character. Full path of the updated mask to be used with the raster file
#' @param NDVI_Thresh numeric. NDVI threshold applied to produce a mask (select pixels with NDVI>NDVI_Thresh)
#' @param Blue_Thresh numeric. Blue threshold applied to produce a mask (select pixels with Blue refl < Blue_Thresh --> filter clouds) refl expected between 0 and 10000
#' @param NIR_Thresh numeric. NIR threshold applied to produce a mask (select pixels with NIR refl < NIR_Thresh) refl expected between 0 and 10000
#' @param Blue numeric. spectral band corresponding to the blue channel (in nanometers)
#' @param Red numeric. spectral band corresponding to the red channel (in nanometers)
#' @param NIR numeric. spectral band corresponding to the NIR channel (in nanometers)
#
# @return MaskPath path for the updated shademask produced
create_mask_from_threshold <- function(ImPath, MaskPath, MaskPath_Update, NDVI_Thresh, Blue_Thresh, NIR_Thresh,
                                       Blue = 480, Red = 690, NIR = 835) {
  # define wavelength corresponding to the spectral domains Blue, Red and NIR
  Spectral_Bands <- c(Blue, Red, NIR)
  ImPathHDR <- get_HDR_name(ImPath)
  HDR <- read_ENVI_header(ImPathHDR)
  # distance between expected bands defining red, blue and NIR info and available band from sensor
  Dist2Band <- 25
  # in case micrometers
  if (!is.null(HDR$`wavelength units`)){
    if (max(HDR$wavelength)<100 | HDR$`wavelength units` == "micrometers"){
      Spectral_Bands <- 0.001*Spectral_Bands
      Dist2Band <- 0.001*Dist2Band
    }
  } else if (is.null(HDR$`wavelength units`)){
    message('wavelength units not provided in the header of the image')
    if (max(HDR$wavelength)<100){
      message('assuming wavelengths are expressed in micrometers')
      Spectral_Bands <- 0.001*Spectral_Bands
      Dist2Band <- 0.001*Dist2Band
    } else {
      message('assuming wavelengths are expressed in nanometers')
    }
  }

  # get image bands correponding to spectral bands of interest
  Image_Bands <- get_image_bands(Spectral_Bands, HDR$wavelength)
  # read band data from image
  Image_Subset <- read_image_bands(ImPath = ImPath, HDR = HDR,
                                   ImBand = Image_Bands$ImBand)
  # create mask
  # check if spectral bands required for NDVI exist

  if (Image_Bands$Distance2WL[2] < Dist2Band & Image_Bands$Distance2WL[3] < Dist2Band) {
    NDVI <- ((Image_Subset[, , 3]) - (Image_Subset[, , 2])) / ((Image_Subset[, , 3]) + (Image_Subset[, , 2]))
  } else {
    NDVI <- matrix(1, nrow = HDR$lines, ncol = HDR$samples)
    message("Could not find the spectral bands required to compute NDVI")
  }
  if (Image_Bands$Distance2WL[1] > Dist2Band) {
    Image_Subset[, , 1] <- Blue_Thresh + 0 * Image_Subset[, , 1]
    message("Could not find a spectral band in the blue domain: will not perform filtering based on blue reflectance")
  }
  if (Image_Bands$Distance2WL[3] > 2*Dist2Band) {
    Image_Subset[, , 3] <- NIR_Thresh + 0 * Image_Subset[, , 3]
    message("Could not find a spectral band in the NIR domain: will not perform filtering based on NIR reflectance")
  }
  Mask <- matrix(0, nrow = HDR$lines, ncol = HDR$samples)
  SelPixels <- which(NDVI > NDVI_Thresh & Image_Subset[, , 1] < Blue_Thresh & Image_Subset[, , 3] > NIR_Thresh)
  Mask[SelPixels] <- 1
  # update initial shade mask
  MaskPath <- update_shademask(MaskPath, HDR, Mask, MaskPath_Update)
  list2Remove <- ls()
  rm(list=list2Remove[-which(list2Remove=='MaskPath')])
  gc()
  return(MaskPath)
}
