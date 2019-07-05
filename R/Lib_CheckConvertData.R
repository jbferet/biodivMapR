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
#' @param Raster.Path character. Full path for the raster to be converted
#' @param Sensor character. Name of the sensor
#' @param Convert.Integer boolean. Should data be converted into integer ?
#' @param Multiplying.Factor numeric. Multiplying factor (eg convert real reflectance values between 0 and 1 into integer between 0 and 10000).
#' @param Output.Directory character. Path to output directory.
#' @param Multiplying.Factor.Last numeric. Multiplying factor for last band.
#'
#' @return Output.Path path for the image converted into ENVI BIL format
#' @import raster
#' @import tools
#' @export
raster2BIL <- function(Raster.Path, Sensor = "unknown", Output.Directory = FALSE, Convert.Integer = TRUE, Multiplying.Factor = 1, Multiplying.Factor.Last = 1) {

  # get directory and file name of original image
  Input.File <- basename(Raster.Path)
  Input.Dir <- dirname(Raster.Path)
  # define path where data will be stored
  if (Output.Directory == FALSE) {
    Output.Path <- paste(Input.Dir, "/Converted_DivMAP/", file_path_sans_ext(Input.File), sep = "")
  } else {
    dir.create(Output.Directory, showWarnings = FALSE, recursive = TRUE)
    Output.Path <- paste(Output.Directory, "/", file_path_sans_ext(Input.File), sep = "")
  }
  message("The converted file will be written in the following location:")
  print(Output.Path)
  Output.Dir <- dirname(Output.Path)
  dir.create(Output.Dir, showWarnings = FALSE, recursive = TRUE)

  # apply multiplying factors
  message("reading initial file")
  Output.Img <- Multiplying.Factor * brick(Raster.Path)
  Last.Band.Name <- Output.Img@data@names[length(Output.Img@data@names)]
  Output.Img[[Last.Band.Name]] <- Multiplying.Factor.Last * Output.Img[[Last.Band.Name]]

  # convert into integer
  if (Convert.Integer == TRUE) {
    Output.Img <- round(Output.Img)
  }

  # write raster
  message("writing converted file")
  if (Convert.Integer == TRUE) {
    r <- writeRaster(Output.Img, filename = Output.Path, format = "EHdr", overwrite = TRUE, datatype = "INT2S")
  } else {
    r <- writeRaster(Output.Img, filename = Output.Path, format = "EHdr", overwrite = TRUE)
  }
  hdr(r, format = "ENVI")

  # remove unnecessary files
  File2Remove <- paste(Output.Path, ".aux.xml", sep = "")
  File2Remove2 <- paste(file_path_sans_ext(Output.Path), ".aux.xml", sep = "")
  file.remove(File2Remove)

  File2Remove <- paste(Output.Path, ".prj", sep = "")
  File2Remove2 <- paste(file_path_sans_ext(Output.Path), ".prj", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Remove <- paste(Output.Path, ".sta", sep = "")
  File2Remove <- paste(file_path_sans_ext(Output.Path), ".sta", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Remove <- paste(Output.Path, ".stx", sep = "")
  File2Remove2 <- paste(file_path_sans_ext(Output.Path), ".stx", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Rename <- paste(file_path_sans_ext(Output.Path), ".hdr", sep = "")
  File2Rename2 <- paste(Output.Path, ".hdr", sep = "")
  if (file.exists(File2Rename)) {
    file.rename(from = File2Rename, to = File2Rename2)
  }

  # change dot into underscore
  Output.Path.US <- file.path(
    dirname(Output.Path),
    gsub(basename(Output.Path), pattern = "[.]", replacement = "_")
  )
  if (!Output.Path.US == Output.Path) {
    file.rename(from = Output.Path, to = Output.Path.US)
  }

  Output.Path.US.HDR <- paste0(Output.Path.US, ".hdr")
  if (!Output.Path.US.HDR == paste0(Output.Path, ".hdr")) {
    file.rename(from = paste0(Output.Path, ".hdr"), to = Output.Path.US.HDR)
    ### UTILITY?? ###
    file.rename(from = Output.Path.US.HDR, to = Output.Path.US.HDR)
  }


  if (!Sensor == "unknown") {
    HDR.Temp.Path <- system.file("extdata", "HDR", paste0(Sensor, ".hdr"), package = "biodivMapR")
    if (file.exists(HDR.Temp.Path)) {
      message("reading header template corresponding to the sensor located here:")
      print(HDR.Temp.Path)
      # get raster template corresponding to the sensor
      HDR.Template <- read.ENVI.header(HDR.Temp.Path)
      # get info to update hdr file
      # read hdr
      HDR.input <- read.ENVI.header(Get.HDR.Name(Output.Path))
      if (!is.null(HDR.Template$wavelength)) {
        HDR.input$wavelength <- HDR.Template$wavelength
      }
      if (!is.null(HDR.Template$`sensor type`)) {
        HDR.input$`sensor type` <- HDR.Template$`sensor type`
      }
      if (!is.null(HDR.Template$`band names`)) {
        HDR.input$`band names` <- HDR.Template$`band names`
      }
      if (!is.null(HDR.Template$`wavelength units`)) {
        HDR.input$`wavelength units` <- HDR.Template$`wavelength units`
      }
      # define visual stretch in the VIS domain
      HDR.input$`default stretch` <- "0 1000 linear"
      # write corresponding hdr file
      write.ENVI.header(HDR.input, Get.HDR.Name(Output.Path))
    } else if (!file.exists(HDR.Temp.Path)) {
      message("Header template corresponding to the sensor expected to be found here")
      print(HDR.Temp.Path)
      message("please provide this header template in order to write info in HDR file")
      print(Get.HDR.Name(Output.Path))
      message("or manually add wavelength location in HDR file, if relevant")
    }
  } else if (Sensor == "unknown") {
    message("please make sure that the follozing header file contains information required")
    print(Get.HDR.Name(Output.Path))
    message("or manually add wavelength location in HDR file, if relevant")
  }
  return(Output.Path)
}

#' Checks if the data to be processed has the format type expected
#'
#' @param Raster.Path full path for the raster to be converted
#' @param Mask is the raster a mask?
#'
#' @export
check_data <- function(Raster.Path, Mask = FALSE) {
  HDR.Path <- Get.HDR.Name(Raster.Path)
  # check if the hdr file exists
  if (file.exists(HDR.Path)) {
    HDR <- read.ENVI.header(HDR.Path)
    if (Mask == FALSE & (!HDR$interleave == "bil") & (!HDR$interleave == "BIL")) {
      message("")
      message("*********************************************************")
      message("The image format may not compatible with the processing chain")
      message("Image format expected:")
      message("ENVI hdr file with band interleaved by line (BIL) file format")
      message("")
      message("Current Image format")
      print(HDR$interleave)
      message("Please run the function named ")
      print("Convert.Raster2BIL")
      message("in order to convert your raster data")
      message("or use appropriate software")
      message("*********************************************************")
      message("")
      stop()
    } else if (Mask == FALSE & ((HDR$interleave == "bil") | (HDR$interleave == "BIL"))) {
      if (HDR$`wavelength units` == "Unknown") {
        message("")
        message("*********************************************************")
        message("Please make sure the wavelengths are in nanometers")
        message("if not, stop processing and convert wavelengths in nanometers in HDR file")
        message("*********************************************************")
        message("")
      }
      if ((!HDR$`wavelength units` == "Nanometers") & (!HDR$`wavelength units` == "nanometers")) {
        message("")
        message("*********************************************************")
        message("Please make sure the wavelengths are in nanometers")
        message("if not, stop processing and convert wavelengths in nanometers in HDR file")
        message("*********************************************************")
        message("")
      }
      if (HDR$`wavelength units` == "micrometers") {
        message("")
        message("*********************************************************")
        message("Please convert wavelengths in nanometers in HDR file")
        message("*********************************************************")
        message("")
        stop()
      }
      if ((HDR$`wavelength units` == "nanometers") | (HDR$`wavelength units` == "Nanometers")) {
        message("")
        message("*********************************************************")
        message("       All information seem OK for image processing      ")
        message("*********************************************************")
        message("")
      }
    }
  } else {
    message("")
    message("*********************************************************")
    message("The following HDR file was expected, but could not be found:")
    print(HDR.Path)
    message("The image format may not compatible with the processing chain")
    message("Image format expected:")
    message("ENVI hdr file with band interleaved by line (BIL) file format")
    message("")
    message("Please run the function named ")
    print("Convert.Raster2BIL")
    message("in order to convert your raster data")
    message("*********************************************************")
    message("")
    stop()
  }
  return(invisible())
}
