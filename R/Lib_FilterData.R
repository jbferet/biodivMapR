# ==============================================================================
# biodivMapR
# Lib_FilterData.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
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
#' @param NIR_Thresh numeric. NIR threshold applied to produce a mask (select pixels with NIR refl < NIR_Thresh) refl expected between 0 and 10000
#' @param Blue numeric. central wavelength corresponding to the blue spectral band (in nanometers)
#' @param Red numeric. central wavelength corresponding to the red spectral band (in nanometers)
#' @param NIR numeric. central wavelength corresponding to the NIR spectral band (in nanometers)
#'
#' @return MaskPath = updated mask file
#' @export
perform_radiometric_filtering <- function(Image_Path, Mask_Path, Output_Dir, TypePCA = "SPCA", NDVI_Thresh = 0.5, Blue_Thresh = 500, NIR_Thresh = 1500, Blue = 480, Red = 700, NIR = 835) {
  # check if format of raster data is as expected
  check_data(Image_Path)
  if (!Mask_Path==FALSE){
    check_data(Mask_Path)
  }
  # define full output directory
  Output_Dir_Full <- define_output_directory(Output_Dir, Image_Path, TypePCA)
  # define dimensions of the image
  ImPathHDR <- get_HDR_name(Image_Path)
  HDR <- read_ENVI_header(ImPathHDR)
  Image_Format <- ENVI_type2bytes(HDR)
  ipix <- as.double(HDR$lines)
  jpix <- as.double(HDR$samples)
  nbPixels <- ipix * jpix
  lenTot <- nbPixels * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image_Format$Bytes) / (1024^3)

  # Create / Update shade mask if optical data
  if (Mask_Path == FALSE | Mask_Path == "") {
    print("Create mask based on NDVI, NIR and Blue threshold")
  } else {
    print("Update mask based on NDVI, NIR and Blue threshold")
  }
  Shade.Update <- paste(Output_Dir_Full, "ShadeMask_Update", sep = "")
  Mask_Path <- create_mask_from_threshold(Image_Path, Mask_Path, Shade.Update, NDVI_Thresh, Blue_Thresh, NIR_Thresh, Blue, Red, NIR)
  return(Mask_Path)
}

# create a mask based on NDVI, Green reflectance and NIR reflectance
# NDVI (min) threshold eliminates non vegetated pixels
# Blue (max) threshold eliminates Clouds
# NIR (min) threshold eliminates shadows
# ! only valid if Optical data!!
#
# @param ImPath full path of a raster file
# @param MaskPath full path of the raster mask corresponding to the raster file
# @param MaskPath.Update wavelength (nm) of the spectral bands to be found
# @param NDVI_Thresh NDVI threshold applied to produce a mask (select pixels with NDVI>NDVI_Thresh)
# @param Blue_Thresh Blue threshold applied to produce a mask (select pixels with Blue refl < Blue_Thresh --> filter clouds) refl expected between 0 and 10000
# @param NIR_Thresh NIR threshold applied to produce a mask (select pixels with NIR refl < NIR_Thresh) refl expected between 0 and 10000
#
# @return MaskPath path for the updated shademask produced
create_mask_from_threshold <- function(ImPath, MaskPath, MaskPath.Update, NDVI_Thresh, Blue_Thresh, NIR_Thresh, Blue = 480, Red = 700, NIR = 835) {
  # define wavelength corresponding to the spectral domains Blue, Red and NIR
  Spectral_Bands <- c(Blue, Red, NIR)
  ImPathHDR <- get_HDR_name(ImPath)
  Header <- read_ENVI_header(ImPathHDR)
  # get image bands correponding to spectral bands of interest
  Image_Bands <- get_image_bands(Spectral_Bands, Header$wavelength)
  # read band data from image
  Image_Subset <- read_image_bands(ImPath, Header, Image_Bands$ImBand)
  # create mask
  # check if spectral bands required for NDVI exist
  if (Image_Bands$Distance2WL[2] < 25 & Image_Bands$Distance2WL[3] < 25) {
    NDVI <- ((Image_Subset[, , 3]) - (Image_Subset[, , 2])) / ((Image_Subset[, , 3]) + (Image_Subset[, , 2]))
  } else {
    NDVI <- matrix(1, nrow = Header$lines, ncol = Header$samples)
    message("Could not find the spectral bands required to compute NDVI")
  }
  if (Image_Bands$Distance2WL[1] > 25) {
    Image_Subset[, , 1] <- Blue_Thresh + 0 * Image_Subset[, , 1]
    message("Could not find a spectral band in the blue domain: will not perform filtering based on blue reflectance")
  }
  if (Image_Bands$Distance2WL[3] > 50) {
    Image_Subset[, , 3] <- NIR_Thresh + 0 * Image_Subset[, , 3]
    message("Could not find a spectral band in the NIR domain: will not perform filtering based on NIR reflectance")
  }
  Mask <- matrix(0, nrow = Header$lines, ncol = Header$samples)
  SelPixels <- which(NDVI > NDVI_Thresh & Image_Subset[, , 1] < Blue_Thresh & Image_Subset[, , 3] > NIR_Thresh)
  Mask[SelPixels] <- 1
  # update initial shade mask
  MaskPath <- update_shademask(MaskPath, Header, Mask, MaskPath.Update)
  return(MaskPath)
}
