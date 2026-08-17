library(biodivMapR)
test_that("spectral_species_maps", {

  library(terra)
  set.seed(123)
  output_dir <- './'
  path_si <- list('CCCI' = system.file("extdata", "RASTER/CCCI.tif", package="biodivMapR"),
                  'NDVI' = system.file("extdata", "RASTER/NDVI.tif", package="biodivMapR"))
  options <- set_options_biodivMapR(fun = 'spectral_species_full')
  options$nb_samples_alpha <- 10000
  spectral_species_full(input_raster_path = path_si,
                        output_dir = output_dir)

  path_res_expected <- system.file("extdata", "RASTER/spectral_species.tiff", package="biodivMapR")
  path_res <- file.path(output_dir, 'spectral_species.tiff')

  expected <- terra::rast(path_res_expected)
  computed <- terra::rast(path_res)
  testthat::expect_equal(c(terra::values(expected)), c(terra::values(computed)))
})
