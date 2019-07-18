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
#' @param Image.Path Path of the image to be processed
#' @param Mask.Path Path of the mask corresponding to the image
#' @param Output.Dir output directory
#' @param TypePCA Type of PCA: "PCA" or "SPCA"
#' @param NDVI.Thresh NDVI threshold applied to produce a mask (select pixels with NDVI>NDVI.Thresh)
#' @param Blue.Thresh Blue threshold applied to produce a mask (select pixels with Blue refl < Blue.Thresh --> filter clouds) refl expected between 0 and 10000
#' @param NIR.Thresh NIR threshold applied to produce a mask (select pixels with NIR refl < NIR.Thresh) refl expected between 0 and 10000
#' @param Mask.Path
#'
#' @return ImPathShade = updated shademask file
#' @export
perform_radiometric_filtering <- function(Image.Path, Mask.Path, Output.Dir, TypePCA = "SPCA", NDVI.Thresh = 0.5, Blue.Thresh = 500, NIR.Thresh = 1500, Blue = 480, Red = 700, NIR = 835) {
  # define full output directory
  Output.Dir.Full <- define_output_directory(Output.Dir, Image.Path, TypePCA)
  # define dimensions of the image
  ImPathHDR <- get_HDR_name(Image.Path)
  HDR <- read_ENVI_header(ImPathHDR)
  Image.Format <- ENVI_type2bytes(HDR)
  ipix <- as.double(HDR$lines)
  jpix <- as.double(HDR$samples)
  Nb.Pixels <- ipix * jpix
  lenTot <- Nb.Pixels * as.double(HDR$bands)
  ImSizeGb <- (lenTot * Image.Format$Bytes) / (1024^3)

  # Create / Update shade mask if optical data
  if (Mask.Path == FALSE | Mask.Path == "") {
    print("Create mask based on NDVI, NIR and Blue threshold")
  } else {
    print("Update mask based on NDVI, NIR and Blue threshold")
  }
  Shade.Update <- paste(Output.Dir.Full, "ShadeMask_Update", sep = "")
  Mask.Path <- create_mask_from_threshold(Image.Path, Mask.Path, Shade.Update, NDVI.Thresh, Blue.Thresh, NIR.Thresh, Blue, Red, NIR)
  return(Mask.Path)
}

# create a mask based on NDVI, Green reflectance and NIR reflectance
# NDVI (min) threshold eliminates non vegetated pixels
# Blue (max) threshold eliminates Clouds
# NIR (min) threshold eliminates shadows
# ! only valid if Optical data!!
#
# @param ImPath full path of a raster file
# @param ImPathShade full path of the raster mask corresponding to the raster file
# @param ImPathShade.Update wavelength (nm) of the spectral bands to be found
# @param NDVI.Thresh NDVI threshold applied to produce a mask (select pixels with NDVI>NDVI.Thresh)
# @param Blue.Thresh Blue threshold applied to produce a mask (select pixels with Blue refl < Blue.Thresh --> filter clouds) refl expected between 0 and 10000
# @param NIR.Thresh NIR threshold applied to produce a mask (select pixels with NIR refl < NIR.Thresh) refl expected between 0 and 10000
#
# @return ImPathShade path for the updated shademask produced
create_mask_from_threshold <- function(ImPath, ImPathShade, ImPathShade.Update, NDVI.Thresh, Blue.Thresh, NIR.Thresh, Blue = 480, Red = 700, NIR = 835) {
  # define wavelength corresponding to the spectral domains Blue, Red and NIR
  Spectral.Bands <- c(Blue, Red, NIR)
  ImPathHDR <- get_HDR_name(ImPath)
  Header <- read_ENVI_header(ImPathHDR)
  # get image bands correponding to spectral bands of interest
  Image.Bands <- get_image_bands(Spectral.Bands, Header$wavelength)
  # read band data from image
  Image.Subset <- read_image_bands(ImPath, Header, Image.Bands$ImBand)
  # create mask
  # check if spectral bands required for NDVI exist
  if (Image.Bands$Distance2WL[2] < 25 & Image.Bands$Distance2WL[3] < 25) {
    NDVI <- ((Image.Subset[, , 3]) - (Image.Subset[, , 2])) / ((Image.Subset[, , 3]) + (Image.Subset[, , 2]))
  } else {
    NDVI <- matrix(1, nrow = Header$lines, ncol = Header$samples)
    message("Could not find the spectral bands required to compute NDVI")
  }
  if (Image.Bands$Distance2WL[1] > 25) {
    Image.Subset[, , 1] <- Blue.Thresh + 0 * Image.Subset[, , 1]
    message("Could not find a spectral band in the blue domain: will not perform filtering based on blue reflectance")
  }
  if (Image.Bands$Distance2WL[3] > 50) {
    Image.Subset[, , 3] <- NIR.Thresh + 0 * Image.Subset[, , 3]
    message("Could not find a spectral band in the NIR domain: will not perform filtering based on NIR reflectance")
  }
  Mask <- matrix(0, nrow = Header$lines, ncol = Header$samples)
  SelPixels <- which(NDVI > NDVI.Thresh & Image.Subset[, , 1] < Blue.Thresh & Image.Subset[, , 3] > NIR.Thresh)
  Mask[SelPixels] <- 1
  # update initial shade mask
  ImPathShade <- update_shademask(ImPathShade, Header, Mask, ImPathShade.Update)
  return(ImPathShade)
}
