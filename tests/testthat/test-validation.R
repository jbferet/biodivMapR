test_that("validation", {

  library(sf)
  library(terra)
  input_raster_path <- list('CCCI' = system.file("extdata", "RASTER/CCCI.tif", package="biodivMapR"),
                            'NDVI' = system.file("extdata", "RASTER/NDVI.tif", package="biodivMapR"))
  window_size <- 10
  input_raster_path_tmp <- system.file("extdata", "RASTER/spectral_species.tiff", package="biodivMapR")
  alpha_metrics <- c('richness', 'shannon', 'simpson', 'hill')
  beta_metrics <- c('bray', 'brayturn', 'jaccard', 'jaccardturn', 'sorensen', 'simpson_diss')

  output_dir <- './'
  set.seed(123)
  input_rast <- terra::rast(input_raster_path_tmp)
  extent_area <- get_raster_extent(input_rast)
  bufferSize <- 0.5*window_size*terra::res(input_rast)[1]
  extent_area_buff <- terra::buffer(x = extent_area, width = -bufferSize)
  nb_samples <- 20
  samples <- sf::st_sample(x = sf::st_as_sf(extent_area_buff), size = nb_samples,
                           force = TRUE)
  # extent_area_buff <- terra::buffer(x = extent_area, width = bufferSize)

  samples <- terra::buffer(x = terra::vect(samples), width = bufferSize,
                           quadsegs = 8, capstyle  = 'square')

  load(system.file("extdata", "RASTER/Kmeans_info.RData", package="biodivMapR"))
  load(system.file("extdata", "RASTER/Beta_info.RData", package="biodivMapR"))
  div <- get_diversity_from_plots(input_rast = terra::rast(unlist(input_raster_path)),
                                  validation_vect = samples,
                                  Kmeans_info = Kmeans_info,
                                  Beta_info = Beta_info,
                                  alpha_metrics = alpha_metrics,
                                  beta_metrics = beta_metrics)

  div$specdiv$source <- NULL
  val <- round(div$specdiv, digits = 5)
  val_expected <- read.table(file = system.file("extdata", "RASTER/validation.csv", package="biodivMapR"),
                             sep = '\t')
  for (valvar in c('richness_mean', 'shannon_mean', 'bray_plot_pcoa_1'))
    testthat::expect_equal(sum(val[[valvar]]), sum(val_expected[[valvar]]))

})
