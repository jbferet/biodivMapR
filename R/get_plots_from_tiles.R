#' get samples for beta diversity mapping
#'
#' @param plotID list.
#' @param plots2sel numeric.
#' @param listfiles character.
#' @param feat_list list.
#' @param window_size numeric.
#' @param minSun numeric.
#' @param p list.
#'
#' @return samples_beta
#' @importFrom dplyr group_split %>%
#' @importFrom terra spatSample extract
#' @importFrom sf st_as_sf st_sf st_crs
#' @export
#'
get_plots_from_tiles <- function(plotID, plots2sel, listfiles, feat_list,
                                 window_size, minSun = 0.75, p = NULL){

  plotID2 <- paste0('_',plotID,'_')
  # tileSI <- listfiles[stringr::str_detect(string = listfiles, pattern = plotID2)]
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
      # read rasters, select more plots than expected to compensate for shaded areas
      selplot <- terra::spatSample(x = rastID, size = as.numeric(plots2sel)+5,
                                   method = "random", na.rm = T, as.df = F,
                                   as.points = T, xy = F, warn = F)

      # define square plots
      buffer <- 10*window_size/2       # = nb pixels x 10 spatial resolution
      circ_plot <- terra::buffer(x = selplot, width = buffer)
      plots_bbox <- plots_square <- list()
      for (i in seq_len(length(selplot))){
        plots_bbox[[i]] <- sf::st_as_sf(circ_plot[i,]) |>
          sf::st_bbox()
        plots_square[[i]] <- sf::st_sf(preprocS2::bbox_to_poly(x = plots_bbox[[i]], crs = sf::st_crs(rastID)))
      }
      df <- do.call(rbind, plots_square)
      vsp <- methods::as(df, "Spatial")
      samples_beta_terra <- terra::vect(vsp)
      # extract info from plots
      samples_beta <- terra::extract(x = rastID, y = samples_beta_terra)
      samples_beta$plotID <- plotID
      samples_beta <- samples_beta %>% group_split(ID)
      # select plots with sufficient sunlit pixels
      nbPixSunLit <- lapply(lapply(samples_beta, '[[', 'mask'), sum, na.rm = T)
      minPix <- minSun*window_size**2
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
  if (!is.null(p)) p()
  return(samples_beta)
}
