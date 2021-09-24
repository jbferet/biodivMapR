test_that("Radiometric Filtering", {

  library(utils)
  destfile <- 'S2A_T33NUD_20180104_Subset'
  # url for the S2 subset
  url <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset'
  download.file(url = url, destfile = destfile, method = 'auto', quiet = FALSE, mode = "w",
                cacheOK = TRUE,
                extra = getOption("download.file.extra"),
                headers = NULL)

  # name your raster HDR
  destfile_HDR <- get_HDR_name(destfile,showWarnings = FALSE)
  # url for the S2 subset header
  urlhdr <-  'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset.hdr'
  download.file(url = urlhdr, destfile = destfile_HDR, method = 'auto', quiet = FALSE, mode = "w",
                cacheOK = TRUE,
                extra = getOption("download.file.extra"),
                headers = NULL)
  Input_Image_File <- destfile
  Input_Mask_File <- FALSE
  Output_Dir <- 'RESULTS'
  window_size <- 10
  nbCPU <- 2
  MaxRAM <- 0.5
  nbclusters <- 50
  NDVI_Thresh <- 0.5
  # these values are relevant only if reflectance is coded as integer values between 0 and 10000.
  Blue_Thresh <- 500
  NIR_Thresh <- 1500
  TypePCA <- 'SPCA'
  print("PERFORM RADIOMETRIC FILTERING")
  Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = Input_Mask_File,
                                                   Output_Dir = Output_Dir, TypePCA = TypePCA,
                                                   NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,NIR_Thresh = NIR_Thresh)
  expect_equal(Input_Mask_File, 'RESULTS/S2A_T33NUD_20180104_Subset/SPCA/ShadeMask_Update')
})
