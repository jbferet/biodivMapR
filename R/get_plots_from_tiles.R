#' get samples for beta diversity mapping
#'
#' @param plotID list.
#' @param plots2sel numeric.
#' @param listfiles character.
#' @param feat_list list.
#' @param window_size numeric.
#' @param min_sun numeric.
#' @param p list.
#'
#' @return samples_beta
#' @importFrom dplyr group_split %>%
#' @importFrom terra spatSample extract
#' @importFrom sf st_as_sf st_sf st_crs
#' @export
#'
get_plots_from_tiles <- function(plotID, plots2sel, listfiles, feat_list,
                                 window_size, min_sun = 0.75, p = NULL){

  plotID2 <- paste0('_',plotID,'_')
  tileSI <- listfiles[grepl(x = listfiles, pattern = plotID2)]
  # get plots for beta diversity
  selpix <- samples_beta <- ID <- NULL
  if (length(tileSI) > 0 & plots2sel >0){
    if (all(file.exists(tileSI))){
      # read rasters
      rastID <- terra::rast(tileSI)
      for (feat in feat_list){
        # whichfeat <- which(stringr::str_detect(basename(terra::sources(rastID)), feat) )
        whichfeat <- which(grepl(x = basename(terra::sources(rastID)),
                                 pattern = feat))
        names(rastID)[whichfeat] <- feat
      }
      # read rasters, select more plots to compensate for shaded areas
      selplot <- terra::spatSample(x = rastID, size = as.numeric(plots2sel)+5,
                                   method = "random", na.rm = TRUE,
                                   as.df = FALSE, as.points = TRUE, xy = FALSE,
                                   warn = FALSE)

      extent_area <- get_raster_extent(rastID[[1]])
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
        raster_res_init <- terra::res(rastID[[1]])[1]
        # raster resolution in meters
        raster_res <- raster_res_init/deg_to_meters
      } else {
        raster_res <- terra::res(rastID[[1]])[1]
      }
      # define square plots
      buff <- raster_res*window_size/2       # = nb pixels x spatial resolution
      circ_plot <- terra::buffer(x = selplot, width = buff)
      plots_bbox <- plots_square <- list()
      for (i in seq_len(length(selplot))){
        plots_bbox[[i]] <- sf::st_as_sf(circ_plot[i,]) |>
          sf::st_bbox()
        plots_square[[i]] <- sf::st_sf(preprocS2::bbox_to_poly(x = plots_bbox[[i]],
                                                               crs = sf::st_crs(rastID)))
      }
      df <- do.call(rbind, plots_square)
      vsp <- methods::as(df, "Spatial")
      samples_beta_terra <- terra::vect(vsp)
      # extract info from plots
      samples_beta <- terra::extract(x = rastID, y = samples_beta_terra)
      samples_beta$plotID <- plotID
      samples_beta <- samples_beta %>% group_split(ID)
      # select plots with sufficient sunlit pixels
      nbPixSunLit <- lapply(lapply(samples_beta, '[[', 'mask'),
                            sum, na.rm = TRUE)
      minPix <- min_sun*window_size**2
      selPlots <- which(unlist(nbPixSunLit)>minPix)
      samples_beta <- samples_beta[selPlots[1:plots2sel]]
      if (length(samples_beta)>0){
        for (i in 1:length(samples_beta)){
          elim <- which(is.na(samples_beta[[i]]$mask))
          if (length(elim)>0)
            samples_beta[[i]] <- samples_beta[[i]][-elim, ]
        }
      } else {
        samples_beta < NULL
      }
    }
  }
  if (!is.null(p))
    p()
  return(samples_beta)
}
