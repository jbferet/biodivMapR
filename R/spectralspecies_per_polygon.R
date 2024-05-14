#' get spectral species corresponding to polygns in a SpatVector object
#
#' @param SpatVector SpatVector.
#' @param input_rast SpatRaster.
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param SelectBands numeric. bands selected from input data
#' @param input_mask SpatRaster.
#' @param MinSun numeric. minimum amount of sunlit pixels in the plots
#
#' @return list
#' @importFrom dplyr group_by
#' @importFrom tidyr nest
#' @importFrom terra extract sources
#' @export
#'
spectralspecies_per_polygon <- function(SpatVector, input_rast,
                                        Kmeans_info, SelectBands = NULL,
                                        input_mask = NULL, MinSun = 0.25){
  # extract pixel info from vector data
  nbPlots_init <- length(SpatVector)
  SpatVector$ID_biodivMapR <- seq_len(nbPlots_init)
  rast_sample <- sample_raster(input_rast = input_rast,
                               pix2extract = SpatVector)
  # get number of pixels per plot
  nbPix_per_plot_init <- data.frame(table(rast_sample$ID))
  # account for mask if provided
  if (!is.null(input_mask)){
    mask_sample <- terra::extract(x = input_mask, y = SpatVector)
    sel <- which(mask_sample[,2]==1)
    rast_sample <- rast_sample[sel,]
  }
  rast_sample <- clean_NAsInf(rast_sample)
  # get plot size
  nbPix_per_plot <- data.frame(table(rast_sample$ID))
  # only get common plots between nbPix_per_plot and nbPix_per_plot nbPix_per_plot_init
  if (length(nbPix_per_plot$Freq)>0){
    whichPlots2keep <- as.numeric(intersect(nbPix_per_plot[,'Var1'],
                                            nbPix_per_plot_init[,'Var1']))

    nbPix_per_plot <- nbPix_per_plot[which(nbPix_per_plot[,'Var1']%in%whichPlots2keep),]
    nbPix_per_plot_init <- nbPix_per_plot_init[which(nbPix_per_plot_init[,'Var1']%in%whichPlots2keep),]
    pcSunlit <- nbPix_per_plot$Freq/nbPix_per_plot_init$Freq
    Sunlit_plots <- which(pcSunlit>=MinSun)
    nbPix_per_plot <- nbPix_per_plot[Sunlit_plots,]
    nbPix_per_plot_init <- nbPix_per_plot_init[Sunlit_plots,]
    nbPlots <- length(nbPix_per_plot_init[,'Var1'])

    selpix <- which(rast_sample$ID%in%nbPix_per_plot_init[,'Var1'])
    rast_sample <- rast_sample[selpix,]
    if (is.null(SelectBands)) SelectBands <- seq_len(dim(rast_sample)[2]-1)
    rast_sample_noID <- rast_sample
    rast_sample_noID$ID <- NULL
    SSValid <- get_spectralSpecies(inputdata = rast_sample_noID,
                                   Kmeans_info = Kmeans_info,
                                   SelectBands = SelectBands)
    SSValid$win_ID <- rast_sample$ID
    AttributeTable <- values(SpatVector)
    AttributeTable$source <- basename(terra::sources(SpatVector))
  } else {
    SSValid <-  NULL
    AttributeTable <- values(SpatVector)
    AttributeTable$source <- basename(sources(SpatVector))
  }
  return(list('SSValid' = SSValid,
              'AttributeTable' = AttributeTable))
}
