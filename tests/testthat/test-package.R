test_that("Radiometric Filtering", {

  Input_Image_File <- system.file('extdata', 'RASTER', 'S2A_T33NUD_20180104_Subset', package = 'biodivMapR')
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
  expect_equal(Input_Mask_File, 'RESULTS/S2A_T33NUD_20180104_Subset/SPCA/ShadeMask_Updat')
})

test_that("Radiometric Filtering", {
  TypePCA <- 'SPCA'
  FilterPCA <- TRUE

  Excluded_WL <- c(0, 400)
  Excluded_WL <- rbind(Excluded_WL, c(895, 1005))
  Excluded_WL <- rbind(Excluded_WL, c(1180, 1480))
  Excluded_WL <- rbind(Excluded_WL, c(1780, 2040))

  Continuum_Removal <- TRUE
  # Perform dimensionality reduction
  print("PERFORM DIMENSIONALITY REDUCTION")
  PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File, Input_Mask_File = Input_Mask_File,
                            Output_Dir = Output_Dir, TypePCA = TypePCA, FilterPCA=FilterPCA,
                            nbCPU = nbCPU, MaxRAM = MaxRAM, Continuum_Removal = Continuum_Removal)

})
