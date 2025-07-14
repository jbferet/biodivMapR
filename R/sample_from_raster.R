#' sample pixels or plots from raster data
#'
#' @param extent_area .
#' @param nb_samples numeric. number of samples to be extracted
#' @param input_rast SpatRaster. raster to extract data from
#' @param input_mask SpatRaster. mask corresponding to raster to extract data from
#' @param window_size numeric. window size for square plots
#' @param capstyle character. shape of the plot, see terra::buffer for info on options for capstyle
#'
#' @return rast_sample dataframe. pixel/plot info extracted from input_rast
#' @importFrom sf st_sample st_as_sf st_crs st_sfc st_crop
#' @importFrom terra vect buffer extract res linearUnits project centroids ext crs
#' @importFrom methods as
#' @importFrom stats runif
#' @importFrom preprocS2 meters_to_decdeg
#' @importFrom crsuggest suggest_crs
#' @export

sample_from_raster <- function(extent_area,
                               nb_samples,
                               input_rast,
                               input_mask = NULL,
                               window_size = NULL,
                               capstyle  = 'square'){

  if (nb_samples < 100000 | !is.null(window_size)){
    if (nb_samples > 100000 & !is.null(window_size)){
      message('extracting large number of plots will take some time')
      message('define window_size = NULL to extract pixels only')
    }
    # randomly sample square cells within the image and mask
    # first deal with global extent, if in this case, to avoid sampling outside of
    # world footprint ...
    if (terra::crs(extent_area, describe= T)$code==4326){
      minmax_abs <- list('xmin' = -179.999,
                         'xmax' = 179.999,
                         'ymin' = -89.999,
                         'ymax' = 89.999)
      minmax_ext <- sf::st_bbox(extent_area)
      if (minmax_ext$xmin < minmax_abs$xmin)
        minmax_ext$xmin <- minmax_abs$xmin
      if (minmax_ext$ymin < minmax_abs$ymin)
        minmax_ext$ymin <- minmax_abs$ymin
      if (minmax_ext$xmax > minmax_abs$xmax)
        minmax_ext$xmax <- minmax_abs$xmax
      if (minmax_ext$ymax > minmax_abs$ymax)
        minmax_ext$ymax <- minmax_abs$ymax
      crs_init <- sf::st_crs(extent_area)
      polyg <- rbind(c(minmax_ext$xmin, minmax_ext$ymin),
                     c(minmax_ext$xmax, minmax_ext$ymin),
                     c(minmax_ext$xmax, minmax_ext$ymax),
                     c(minmax_ext$xmin, minmax_ext$ymax),
                     c(minmax_ext$xmin, minmax_ext$ymin)) |>
        list() |>
        sf::st_polygon() |>
        sf::st_sfc(crs = crs_init)
      extent_area <- sf::st_crop(x = sf::st_as_sf(extent_area),  y = polyg)
      extent_area <- as(object = extent_area, Class = 'SpatVector')
    }
    samples <- sf::st_sample(x = sf::st_as_sf(extent_area), size = nb_samples,
                             force = TRUE)

    # define cell size based on window_size
    if (is.null(window_size)) {
      samples <- terra::vect(samples)
    } else {
      # deal with crs units when not meters
      if (sf::st_crs(sf::st_as_sf(extent_area))$units_gdal=='degree'){
        # get resolution in meters for centroid
        centroid_aoi <- terra::centroids(x = extent_area)
        occs_df <- data.frame('latitude' = terra::ext(centroid_aoi)[3],
                              'longitude' = terra::ext(centroid_aoi)[1],
                              'distance' = 1)
        # how many degrees for one meter?
        distlatlon <- preprocS2::meters_to_decdeg(occs_df = occs_df,
                                                  lat_col = 'latitude',
                                                  lon_col = 'longitude',
                                                  distance = 'distance')
        deg_to_meters <- mean(unlist(distlatlon))
        raster_res_init <- terra::res(input_rast[[1]])[1]
        # raster resolution in meters
        raster_res <- raster_res_init/deg_to_meters
      } else {
        raster_res <- terra::res(input_rast[[1]])[1]
      }
      # apply buffer
      bufferSize <- raster_res*window_size/2
      samples <- terra::buffer(x = terra::vect(samples), width = bufferSize,
                               quadsegs = 8, capstyle  = capstyle)
    }
    rast_sample <- sample_raster(input_rast = input_rast, pix2extract = samples)
    # account for mask if provided
    if (!is.null(input_mask)){
      mask_sample <- terra::extract(x = input_mask, y = samples)
      sel <- which(mask_sample[,2]==1)
      rast_sample <- rast_sample[sel,]
    }
    rast_sample <- clean_NAsInf(rast_sample)
  } else {
    xysamples <- get_xy_samples(input_rast, nb_samples, input_mask = input_mask)
    rast_sample <- sample_raster(input_rast = input_rast,
                                 pix2extract = xysamples,
                                 xy = TRUE)
    rast_sample <- clean_NAsInf(rast_sample)
  }
  return(rast_sample)
}
