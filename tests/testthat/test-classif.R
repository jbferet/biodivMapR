test_that("diversity_maps_from_classif", {

  library(terra)
  set.seed(123)
  output_dir <- './'
  input_raster_path <- system.file("extdata", "RASTER/spectral_species.tiff", package="biodivMapR")
  path_ss <- biodivMapR_classif(input_raster_path = input_raster_path,
                                output_dir = output_dir, window_size = 10,
                                nb_samples_beta = 100)

  path_res_expected <- system.file("extdata", "RASTER/shannon_classif.tiff", package="biodivMapR")
  path_res <- file.path(output_dir, 'shannon_classif.tiff')
  expected <- terra::rast(path_res_expected)
  computed <- terra::rast(path_res)
  testthat::expect_equal(c(terra::values(expected)), c(terra::values(computed)))

  # path_res_expected <- system.file("extdata", "RASTER/pcoa_bray_classif.tiff", package="biodivMapR")
  # path_res <- file.path(output_dir, 'pcoa_bray_classif.tiff')
  # expected <- terra::rast(path_res_expected)
  # computed <- terra::rast(path_res)
  # testthat::expect_equal(c(terra::values(expected$pco1)), c(terra::values(computed$pco1)))

})
