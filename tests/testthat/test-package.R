test_that("Radiometric Filtering", {
  library(utils)
  destfile <- 'S2A_T33NUD_20180104_Subset'
  # url for the S2 subset
  url <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset'
  download.file(url = url, destfile = destfile)

  # name your raster HDR
  destfile_HDR <- get_HDR_name(destfile,showWarnings = FALSE)
  # url for the S2 subset header
  urlhdr <-  paste0(url, '.hdr')
  download.file(url = urlhdr, destfile = destfile_HDR)
  output_dir <- 'RESULTS'
  dir.create(output_dir,showWarnings = F, recursive = T)
  S2_hdr_file <- system.file("extdata", "HDR/SENTINEL_2.hdr", package="biodivMapR")
  S2_hdr <- read_ENVI_header(HDRpath = S2_hdr_file)
  updated_mask_path <- radiometric_filtering(input_raster_path = destfile,
                                             output_dir = output_dir,
                                             input_rast_wl = S2_hdr$wavelength)
  testthat::expect_equal(updated_mask_path, 'RESULTS/mask_update.tiff')
})
