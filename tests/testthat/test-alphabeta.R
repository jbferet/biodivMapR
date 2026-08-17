test_that("diversity_maps", {
  library(terra)
  input_raster_path <- list('CCCI' = system.file("extdata", "RASTER/CCCI.tif", package="biodivMapR"),
                            'NDVI' = system.file("extdata", "RASTER/NDVI.tif", package="biodivMapR"))
  set.seed(123)
  window_size <- 10
  alpha_metrics <- c('richness', 'shannon', 'simpson', 'hill')
  beta_metrics <- c('bray', 'brayturn', 'jaccard', 'jaccardturn', 'sorensen', 'simpson_diss')
  # get diversity maps
  output_dir <- './'
  options <- set_options_biodivMapR(fun = 'biodivMapR')
  options$alpha_metrics <- alpha_metrics
  options$beta_metrics <- beta_metrics
  options$nb_samples_beta <- 40
  options$moving_window <- FALSE
  bdm_test<- biodivMapR(input_raster_path = input_raster_path,
                        output_dir = output_dir, window_size = window_size,
                        nbCPU = 1, options = options)


  path_res_expected <- list('shannon' = system.file("extdata", "RASTER/shannon_mean.tiff", package="biodivMapR"),
                            'richness' = system.file("extdata", "RASTER/richness_mean.tiff", package="biodivMapR"),
                            'bray' = system.file("extdata", "RASTER/bray.tiff", package="biodivMapR"),
                            'brayturn' = system.file("extdata", "RASTER/brayturn.tiff", package="biodivMapR"))

  path_res <- list('shannon' = file.path(output_dir, 'shannon_mean.tiff'),
                   'richness' = file.path(output_dir, 'richness_mean.tiff'),
                   'bray' = file.path(output_dir, 'bray.tiff'),
                   'brayturn' = file.path(output_dir, 'brayturn.tiff'))

  corr <- list()
  for (idx in names(path_res)){
    expected <- terra::rast(path_res_expected[[idx]])
    computed <- terra::rast(path_res[[idx]])
    testthat::expect_equal(c(terra::values(expected)), c(terra::values(computed)))
  }
})
