#' extract pixel information from a raster based on vector footprint
#'
#' @param SpatVector SpatVector object SpatRaster or list of SpatRaster
#' @param input_rast SpatRaster input raster
#' @param input_mask SpatRaster optioal mask
#' @param MinSun numeric. minimum sunlit/unmasked proportion
#' @param prog boolean progressbar ?
#'
#' @return list including rast_sample_vect: dataframe corresponding to
#' information extracted from input_rast, and AttributeTable from SpatVector

#' @importFrom terra sources extract
#' @export

extract_vect_from_rast <- function(SpatVector, input_rast,
                                   input_mask = NULL, MinSun = 0.25,
                                   prog = T){
  # extract pixel info from vector data
  SpatVector$ID_biodivMapR <- seq_len(length(SpatVector))
  AttributeTable <- values(SpatVector)
  AttributeTable$source <- basename(terra::sources(SpatVector))
  rast_sample_vect <- sample_raster(input_rast = input_rast,
                                    pix2extract = SpatVector,
                                    prog = prog)
  # get number of pixels per plot from initial footprints
  nbPix_per_plot_init <- data.frame(table(rast_sample_vect$ID))
  # eliminate masked pixels if relevant
  if (!is.null(input_mask)){
    mask_sample <- terra::extract(x = input_mask, y = SpatVector)
    sel <- which(mask_sample[,2]==1)
    rast_sample_vect <- rast_sample_vect[sel,]
  }
  rast_sample_vect <- clean_NAsInf(rast_sample_vect)
  # get plot size
  nbPix_per_plot <- data.frame(table(rast_sample_vect$ID))
  # only get common plots between nbPix_per_plot and nbPix_per_plot nbPix_per_plot_init
  if (length(nbPix_per_plot$Freq)>0){
    whichPlots2keep <- as.numeric(intersect(nbPix_per_plot[,'Var1'],
                                            nbPix_per_plot_init[,'Var1']))

    nbPix_per_plot <- nbPix_per_plot[which(nbPix_per_plot[,'Var1']%in%whichPlots2keep),]
    nbPix_per_plot_init <- nbPix_per_plot_init[which(nbPix_per_plot_init[,'Var1']%in%whichPlots2keep),]
    Sunlit_plots <- which(nbPix_per_plot$Freq/nbPix_per_plot_init$Freq >= MinSun)
    nbPix_per_plot <- nbPix_per_plot[Sunlit_plots,]
    nbPix_per_plot_init <- nbPix_per_plot_init[Sunlit_plots,]
    selpix <- which(rast_sample_vect$ID%in%nbPix_per_plot_init[,'Var1'])
    rast_sample_vect <- rast_sample_vect[selpix,]
  }
  return(list('rast_sample_vect' = rast_sample_vect,
              'AttributeTable' = AttributeTable))
}
