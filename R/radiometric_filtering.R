#' Performs radiometric filtering based on three criteria: NDVI, NIR reflectance, Blue reflectance
#'
#' @param input_raster_path character. Path of the image to be processed
#' @param output_dir character. Path for output directory
#' @param input_mask_path character. Path of the mask corresponding to the image
#' @param input_rast_wl numeric. spectral bands used to identify relevant bands
#' @param NDVI_Thresh numeric. NDVI threshold applied to produce a mask (select pixels with NDVI>NDVI_Thresh)
#' @param Blue_Thresh numeric. Blue threshold applied to produce a mask (select pixels with Blue refl < Blue_Thresh --> filter clouds) refl expected between 0 and 10000
#' @param NIR_Thresh numeric. NIR threshold applied to produce a mask (select pixels with NIR refl > NIR_Thresh) refl expected between 0 and 10000
#' @param Blue numeric. central wavelength corresponding to the blue spectral band (in nanometers)
#' @param Red numeric. central wavelength corresponding to the red spectral band (in nanometers)
#' @param NIR numeric. central wavelength corresponding to the NIR spectral band (in nanometers)
#' @param maxRows numeric. max number of rows in each block
#' @param filetype character. GDAL driver
#' @param maskfilename character. name of the updated mask file
#'
#' @return MaskPath = updated mask file
#' @importFrom terra rast blocks readStart writeStart writeValues readStop
#' @export

radiometric_filtering <- function(input_raster_path, output_dir, input_rast_wl,
                                  input_mask_path = NULL, NDVI_Thresh = 0.65,
                                  Blue_Thresh = 500, NIR_Thresh = 1500,
                                  Blue = 480, Red = 670, NIR = 835,
                                  maxRows = 1000, filetype = 'GTiff',
                                  maskfilename = 'mask_update') {

  # produce SpatRaster
  input_rast <- terra::rast(input_raster_path)
  names(input_rast) <- input_rast_wl
  # check if format of raster data is as expected
  check_data(input_data = input_rast, arguments = 'input_rast')
  if (!is.null(input_mask_path)) {
    input_mask <- terra::rast(input_mask_path)
    check_data(input_data = input_mask, arguments = 'input_mask')
  }
  mask_update <- file.path(output_dir, maskfilename)
  dir.create(path = output_dir,recursive = T,showWarnings = F)
  if (filetype%in%c('GTiff', 'COG') & ! grepl(x = mask_update, pattern = '.tif')) {
    mask_update <- paste0(mask_update, '.tif')
  }

  # wavelengths expected to perform filtering
  Spectral_Bands <- data.frame('Blue' = Blue, 'Red' = Red, 'NIR' = NIR)
  wl_filter <- list('filter' = 'WL', 'Spectral_Bands' = Spectral_Bands)
  check_data(input_data = names(input_rast), arguments = wl_filter)
  # get image bands correponding to spectral bands of interest
  Spectral_Bands <- data.frame('Blue' = Blue,
                               'Red' = Red,
                               'NIR' = NIR)
  Thresholds <- data.frame('Blue' = Blue_Thresh,
                           'NDVI' = NDVI_Thresh,
                           'NIR' = NIR_Thresh)
  Image_Bands <- get_image_bands(Spectral_Bands = Spectral_Bands,
                                 wavelength = as.numeric(names(input_rast)))

  # read image chunks in order to filter data
  r_in <- list()
  for (fid in names(Image_Bands$ImBand)){
    if (!is.na(Image_Bands$ImBand[[fid]])){
      r_in[[fid]] <- input_rast[[Image_Bands$ImBand[[fid]]]]
    }
  }
  if (!is.null(input_mask_path)) r_in$mask <- input_mask
  # Adjust size of individual chunks
  brast <- terra::rast(input_rast[[1]])
  blk <- terra::blocks(brast[[1]], n = 10)
  blk <- maxRows_chunk(blk = blk, maxRows = maxRows)
  for (fid in names(r_in)) terra::readStart(r_in[[fid]])
  # write mask
  r_out <- r_in[[1]]
  names(r_out) <- 'mask_update'
  r2 <- terra::writeStart(x = r_out, filename = mask_update,
                          overwrite = TRUE, filetype = filetype,
                          datatype = "INT1U")

  # produce a list of blocks to read and process
  blk_list <- list()
  for (chunk in seq_len(length(blk$row))) blk_list[[chunk]] <- list('row'= blk$row[chunk],
                                                                  'nrows' = blk$nrows[chunk])
  # compute diversity metrics for each block
  for (bloc in blk_list){
    update_mask <- radiometricfilter_chunk(blk = bloc, r_in = r_in,
                                           Thresholds = Thresholds,
                                           Spectral_Bands = Spectral_Bands)
    terra::writeValues(x = r_out, v = update_mask,
                       bloc$row, nrows = bloc$nrows)
  }
  for (fid in names(r_in)) terra::readStop(r_in[[fid]])
  terra::writeStop(r_out)
  return(mask_update)
}

