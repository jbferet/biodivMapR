#' computes spectral indices over a spatial extent from a raster
#'
#' @param aoi spatVector.
#' @param aoi_ID character.
#' @param rastobj spatRaster
#' @param si_list character.
#' @param output_dir character.
#' @param maskobj spatRaster
#' @param sensor_name character.
#' @param radiometric_filter list.
#' @param p object
#' @param overwrite boolean.
#' @param site_name character.
#'
#' @return xxx
#' @importFrom terra vect crop res writeRaster values
#' @importFrom methods as
#' @importFrom sf st_sf
#' @importFrom spinR compute_S2SI_Raster
#' @export

get_si_tiles_from_raster <- function(aoi, aoi_ID, rastobj, si_list, output_dir,
                                     maskobj = NULL, sensor_name = 'sentinel-2',
                                     radiometric_filter = list('cloudMask' = NULL,
                                                               'shadeMask' = NULL,
                                                               'NDVIMask' = NULL),
                                     p = NULL, overwrite = FALSE, site_name = NULL){

  # define default site_name
  if (is.null(site_name))
    site_name <- 'plot'
  # make sure directories are already created
  output_dir_mask <- file.path(output_dir, 'mask')
  dir.create(path = output_dir_mask, showWarnings = FALSE, recursive = TRUE)
  output_dir_si <- file.path(output_dir, 'spectral_indices')
  dir.create(path = output_dir_si, showWarnings = FALSE, recursive = TRUE)

  # which files are expected as outputs
  filename_mask <- file.path(output_dir_mask, paste0('mask_',aoi_ID,'.tiff'))
  filename_si <- as.list(file.path(output_dir_si,
                                   paste0(site_name, '_', aoi_ID,
                                          '_', si_list, '.tiff')))
  names(filename_si) <- si_list

  # if overwrite or one of expected files does not exists
  if (! file.exists(filename_mask) | ! all(file.exists(unlist(filename_si))) | overwrite){
    aoi_plot <- methods::as(sf::st_sf(aoi), "Spatial")
    aoi_plot <- terra::vect(aoi_plot)
    sensor_refl <- terra::crop(x = rastobj, y = aoi_plot)
    # handles mask
    if (! is.null(maskobj)) {
      bin_mask <- terra::crop(x = maskobj, y = aoi_plot)
      bin_mask <- terra::resample(x = bin_mask, y = sensor_refl,
                                  method  = 'near')
    } else {
      bin_mask <- 0*sensor_refl[[1]]+1
    }
    if ('S2' %in% toupper(sensor_name) | 'SENTINEL2' %in% toupper(sensor_name) |
        'SENTINEL_2' %in% toupper(sensor_name) |
        'SENTINEL-2' %in% toupper(sensor_name)){
      HDRpath <- system.file('extdata', 'HDR', 'Sentinel_2.hdr',
                             package = 'biodivMapR')
      sensor <- 'S2'
    } else if ('LANDSAT' %in% toupper(sensor_name) |
               'LANDSAT7' %in% toupper(sensor_name) |
               'LANDSAT_7' %in% toupper(sensor_name) |
               'LANDSAT-7' %in% toupper(sensor_name) |
               'LANDSAT8' %in% toupper(sensor_name) |
               'LANDSAT_8' %in% toupper(sensor_name) |
               'LANDSAT-8' %in% toupper(sensor_name) |
               'LANDSAT9' %in% toupper(sensor_name) |
               'LANDSAT_9' %in% toupper(sensor_name) |
               'LANDSAT-9' %in% toupper(sensor_name)){
      HDRpath <- system.file('extdata', 'HDR', 'Landsat_7.hdr', package = 'biodivMapR')
      sensor <- 'landsat'
    }
    hdr <- biodivMapR::read_ENVI_header(HDRpath = HDRpath)
    SensorBands <- hdr$wavelength
    if (sensor == 'landsat'){
      bandnames <- strsplit(x = hdr$`band names`, split = ', ')[[1]]
      names(sensor_refl) <- bandnames
      sensor_refl$blue <- sensor_refl$B01
      sensor_refl$red <- sensor_refl$B03
      sensor_refl$nir <- sensor_refl$B04
    } else if (sensor == 'S2'){
      if (terra::res(sensor_refl)[1] == 10)
        bandnames <- c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B08', 'B8A', 'B11', 'B12')
      if (terra::res(sensor_refl)[1] == 20)
        bandnames <- c('B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B8A', 'B11', 'B12')
      names(sensor_refl) <- bandnames
      sensor_refl$blue <- sensor_refl$B02
      sensor_refl$red <- sensor_refl$B04
      if (terra::res(sensor_refl)[1] == 10){
        sensor_refl$nir <- sensor_refl$B08
      } else if (terra::res(sensor_refl)[1] == 20){
        sensor_refl$nir <- sensor_refl$B8A
      }
    }

    # compute radiometric mask
    # set values out of range if NA to keep all pixels
    ndvi <- (sensor_refl$nir-sensor_refl$red)/(sensor_refl$nir+sensor_refl$red)
    if (is.null(radiometric_filter$NDVIMask))
      radiometric_filter$NDVIMask <- 0
    if (is.null(radiometric_filter$cloudMask)){
      sensor_refl$blue <- 0*sensor_refl[[1]]
      radiometric_filter$cloudMask <- 1
    }
    if (is.null(radiometric_filter$shadeMask)){
      sensor_refl$nir   <- 1 + 0*sensor_refl[[1]]
      radiometric_filter$shadeMask <- 0
    }
    sel <- sensor_refl$blue < radiometric_filter$cloudMask &
      sensor_refl$nir > radiometric_filter$shadeMask &
      ndvi > radiometric_filter$NDVIMask
    # mainmask <- preprocS2::get_mainmask(mask_path = , sensor_refl, aoi_plot_sf)
    bin_mask[is.na(bin_mask)] <- 0
    names(bin_mask) <- 'mask'
    bin_mask <- bin_mask*sel
    # save mask
    terra::writeRaster(x = bin_mask, filename = filename_mask, overwrite = TRUE)

    # compute SI
    SI <- spinR::compute_S2SI_Raster(Refl = sensor_refl,
                                     SensorBands = SensorBands,
                                     Sel_Indices = si_list,
                                     StackOut = F,
                                     ReflFactor = 10000)
    # mask SI
    na_mask <- bin_mask
    na_mask[which(terra::values(bin_mask)==0)] <- NA
    # save SI
    for (si in si_list){
      sisel <- SI$SpectralIndices[[si]]
      names(sisel) <- si
      sisel <- sisel*na_mask
      terra::writeRaster(x = sisel,
                         filename = filename_si[[si]],
                         overwrite = TRUE)
    }
  }
  if (!is.null(p))
    p()
  return(output_dir_si)
}

