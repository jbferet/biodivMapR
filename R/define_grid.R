#' define a grid over a raster
#'
#' @param raster_path character. path for raster
#' @param cellsize numeric. size of grid cell
#'
#' @return list
#' @export

define_grid <- function(raster_path, cellsize){
  rast_obj <- terra::rast(raster_path)[[1]]
  mask_polygon <- terra::as.polygons(x = terra::ext(rast_obj),
                                     crs=terra::crs(rast_obj))
  aoi <- sf::st_as_sf(mask_polygon, quiet = T)
  crs_init <- sf::st_crs(aoi)
  aoi <- sf::st_multipolygon(aoi$geometry)
  aoi <- sf::st_sf(sf::st_sfc(aoi))
  sf::st_crs(aoi) <- crs_init
  aoi_grid <- sf::st_make_grid(aoi, cellsize = cellsize, square = TRUE)
  suppressMessages(sf::sf_use_s2(FALSE))
  intersect <- as.data.frame(sf::st_intersects(x = aoi_grid,
                                               aoi))
  suppressMessages(sf::sf_use_s2(TRUE))
  aoi_grid <- aoi_grid[intersect$row.id]
  sf::st_crs(aoi_grid) <- crs_init
  filename <- file.path(dirname(raster_path), paste0("aoi_tiles_",
                                                     cellsize, "m.gpkg"))
  sf::st_write(obj = aoi_grid, dsn = filename, overwrite = T,
               append = F, driver = "GPKG", quiet = T)
  nb_tiles <- length(aoi_grid)
  plots <- preprocS2::get_plot_list(dsn = filename, nbdigits = nchar(nb_tiles))
  return(list('plots' = plots,
              'aoi_grid_path' = filename))
}
