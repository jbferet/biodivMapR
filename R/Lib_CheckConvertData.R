# ==============================================================================
# biodivMapR
# Lib_CheckConvertData.R
# ==============================================================================
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@irstea.fr>
# Copyright 2018/07 Jean-Baptiste FERET
# ==============================================================================
# This library checks if a raster in format expected for diversity mapping
# if not, can convert the raster into BIL format expected for diversity mapping
# ==============================================================================

#' converts a raster into BIL format as expected by DivMapping codes
#'
#' @param Raster_Path character. Full path for the raster to be converted
#' @param Sensor character. Name of the sensor. a .hdr template for the sensor should be provided in extdata/HDR
#' @param Convert_Integer boolean. Should data be converted into integer ?
#' @param Multiplying_Factor numeric. Multiplying factor (eg convert real reflectance values between 0 and 1 into integer between 0 and 10000).
#' @param Output_Dir character. Path to output directory.
#' @param Multiplying_Factor_Last numeric. Multiplying factor for last band.
#'
#' @return Output_Path path for the image converted into ENVI BIL format
#' @import raster
#' @import tools
#' @export
raster2BIL <- function(Raster_Path, Sensor = "unknown", Output_Dir = FALSE, Convert_Integer = TRUE, Multiplying_Factor = 1, Multiplying_Factor_Last = 1) {

  # get directory and file name of original image
  Input_File <- basename(Raster_Path)
  Input_Dir <- dirname(Raster_Path)
  # define path where data will be stored
  if (Output_Dir == FALSE) {
    Output_Path <- paste(Input_Dir, "/Converted_DivMAP/", file_path_sans_ext(Input_File), sep = "")
  } else {
    dir.create(Output_Dir, showWarnings = FALSE, recursive = TRUE)
    Output_Path <- paste(Output_Dir, "/", file_path_sans_ext(Input_File), sep = "")
  }
  message("The converted file will be written in the following location:")
  print(Output_Path)
  Output_Dir <- dirname(Output_Path)
  dir.create(Output_Dir, showWarnings = FALSE, recursive = TRUE)

  # apply multiplying factors
  message("reading initial file")
  Output_Img <- Multiplying_Factor * brick(Raster_Path)
  Last_Band_Name <- Output_Img@data@names[length(Output_Img@data@names)]
  Output_Img[[Last_Band_Name]] <- Multiplying_Factor_Last * Output_Img[[Last_Band_Name]]

  # convert into integer
  if (Convert_Integer == TRUE) {
    Output_Img <- round(Output_Img)
  }

  # write raster
  message("writing converted file")
  if (Convert_Integer == TRUE) {
    r <- writeRaster(Output_Img, filename = Output_Path, format = "EHdr", overwrite = TRUE, datatype = "INT2S")
  } else {
    r <- writeRaster(Output_Img, filename = Output_Path, format = "EHdr", overwrite = TRUE)
  }
  hdr(r, format = "ENVI")

  # remove unnecessary files
  File2Remove <- paste(Output_Path, ".aux.xml", sep = "")
  File2Remove2 <- paste(file_path_sans_ext(Output_Path), ".aux.xml", sep = "")
  file.remove(File2Remove)

  File2Remove <- paste(Output_Path, ".prj", sep = "")
  File2Remove2 <- paste(file_path_sans_ext(Output_Path), ".prj", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Remove <- paste(Output_Path, ".sta", sep = "")
  File2Remove <- paste(file_path_sans_ext(Output_Path), ".sta", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Remove <- paste(Output_Path, ".stx", sep = "")
  File2Remove2 <- paste(file_path_sans_ext(Output_Path), ".stx", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Rename <- paste(file_path_sans_ext(Output_Path), ".hdr", sep = "")
  File2Rename2 <- paste(Output_Path, ".hdr", sep = "")
  if (file.exists(File2Rename)) {
    file.rename(from = File2Rename, to = File2Rename2)
  }

  # change dot into underscore
  Output_Path_US <- file.path(
    dirname(Output_Path),
    gsub(basename(Output_Path), pattern = "[.]", replacement = "_")
  )
  if (!Output_Path_US == Output_Path) {
    file.rename(from = Output_Path, to = Output_Path_US)
  }

  Output_Path_US_HDR <- paste0(Output_Path_US, ".hdr")
  if (!Output_Path_US_HDR == paste0(Output_Path, ".hdr")) {
    file.rename(from = paste0(Output_Path, ".hdr"), to = Output_Path_US_HDR)
    ### UTILITY?? ###
    file.rename(from = Output_Path_US_HDR, to = Output_Path_US_HDR)
  }

  if (!Sensor == "unknown") {
    HDR_Temp_Path <- system.file("extdata", "HDR", paste0(Sensor, ".hdr"), package = "biodivMapR")
    if (file.exists(HDR_Temp_Path)) {
      message("reading header template corresponding to the sensor located here:")
      print(HDR_Temp_Path)
      # get raster template corresponding to the sensor
      HDR_Template <- read_ENVI_header(HDR_Temp_Path)
      # get info to update hdr file
      # read hdr
      HDR_input <- read_ENVI_header(get_HDR_name(Output_Path))
      if (!is.null(HDR_Template$wavelength)) {
        HDR_input$wavelength <- HDR_Template$wavelength
      }
      if (!is.null(HDR_Template$`sensor type`)) {
        HDR_input$`sensor type` <- HDR_Template$`sensor type`
      }
      if (!is.null(HDR_Template$`band names`)) {
        HDR_input$`band names` <- HDR_Template$`band names`
      }
      if (!is.null(HDR_Template$`wavelength units`)) {
        HDR_input$`wavelength units` <- HDR_Template$`wavelength units`
      }
      # define visual stretch in the VIS domain
      HDR_input$`default stretch` <- "0 1000 linear"
      # write corresponding hdr file
      write_ENVI_header(HDR_input, get_HDR_name(Output_Path))
    } else if (!file.exists(HDR_Temp_Path)) {
      message("Header template corresponding to the sensor expected to be found here")
      print(HDR_Temp_Path)
      message("please provide this header template in order to write info in HDR file")
      print(get_HDR_name(Output_Path))
      message("or manually add wavelength location in HDR file, if relevant")
    }
  } else if (Sensor == "unknown") {
    message("please make sure that the follozing header file contains information required")
    print(get_HDR_name(Output_Path))
    message("or manually add wavelength location in HDR file, if relevant")
  }
  return(Output_Path)
}

#' Checks if the data to be processed has the format type expected
#'
#' @param Raster_Path character. full path for the raster to be converted
#' @param Mask boolean. Set true if the raster is a mask
#'
#' @return nothing
#' @export
check_data <- function(Raster_Path, Mask = FALSE) {
  HDR_Path <- get_HDR_name(Raster_Path)
  # check if the hdr file exists
  if (file.exists(HDR_Path)) {
    HDR <- read_ENVI_header(HDR_Path)
    if (Mask == FALSE & (!HDR$interleave == "bil") & (!HDR$interleave == "BIL")) {
      message("*********************************************************")
      message("The image format may not compatible with the processing chain")
      message("Image format expected:")
      message("ENVI hdr file with band interleaved by line (BIL) file format")
      message("")
      message("Current Image format")
      print(HDR$interleave)
      message("Please run the function named ")
      print("raster2BIL")
      message("in order to convert your raster data")
      message("or use appropriate software")
      message("*********************************************************")
      stop()
    } else if (Mask == FALSE & ((HDR$interleave == "bil") | (HDR$interleave == "BIL"))) {
      if (is.null(HDR$`wavelength units`)) {
        message("*********************************************************")
        message("Image to process is not multispectral/hyperspectral image ")
        message("Format is OK, but make sure Continuum_Removal is set to FALSE")
        message("*********************************************************")
      } else {
        if (HDR$`wavelength units` == "Unknown") {
          message("*********************************************************")
          message("IF MULTI / HYPERSPECTRAL DATA: ")
          message("Please make sure the wavelengths are in nanometers")
          message("if not, stop processing and convert wavelengths in nanometers in HDR file")
          message("*********************************************************")
        }
        if ((!HDR$`wavelength units` == "Nanometers") & (!HDR$`wavelength units` == "nanometers")) {
          message("*********************************************************")
          message("IF MULTI / HYPERSPECTRAL DATA: ")
          message("Please make sure the wavelengths are in nanometers")
          message("if not, stop processing and convert wavelengths in nanometers in HDR file")
          message("*********************************************************")
        }
        if (HDR$`wavelength units` == "micrometers") {
          message("*********************************************************")
          message("Please convert wavelengths in nanometers in HDR file")
          message("*********************************************************")
          stop()
        }
        if ((HDR$`wavelength units` == "nanometers") | (HDR$`wavelength units` == "Nanometers")) {
          message("*********************************************************")
          message("      	  Format of main raster OK for processing      	  ")
          message("*********************************************************")
        }
      }
    } else if (Mask == TRUE & HDR$bands == 1 & ((HDR$interleave == "bil") | (HDR$interleave == "BIL") | (HDR$interleave == "bsq") | (HDR$interleave == "BSQ"))) {
      message("*********************************************************")
      message("         Format of mask raster OK for processing         ")
      message("*********************************************************")
	  } else if (Mask == TRUE & HDR$bands > 1) {
      message("*********************************************************")
      message("       Mask raster should contain only one layer         ")
      message("    Please produce a binary mask with a unique layer     ")
      message("*********************************************************")
  		stop()
	  }
  } else {
    message("*********************************************************")
    message("The following HDR file was expected, but could not be found:")
    print(HDR_Path)
    message("The image format may not compatible with the processing chain")
    message("Image format expected:")
    message("ENVI hdr file with band interleaved by line (BIL) file format")
    message("")
    message("Please run the function named ")
    print("raster2BIL")
    message("in order to convert your raster data")
    message("*********************************************************")
    stop()
  }
  return(invisible())
}
