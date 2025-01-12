#' extract pixel information from a raster based on vector footprint
#'
#' @param SpatVector SpatVector. object
#' @param input_rast SpatRaster. input raster
#' @param input_mask SpatRaster. optional mask
#' @param MinSun numeric. minimum sunlit/unmasked proportion
#' @param prog boolean. progressbar ?
#'
#' @return list including rast_sample: dataframe corresponding to
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
  rast_sample <- sample_raster(input_rast = input_rast,
                               pix2extract = SpatVector,
                               prog = prog)
  # get number of pixels per plot from initial footprints
  nbPix_per_plot_init <- data.frame(table(rast_sample$ID))
  # eliminate masked pixels if relevant
  if (!is.null(input_mask)){
    mask_sample <- terra::extract(x = input_mask, y = SpatVector)
    sel <- which(mask_sample[,2]==1)
    rast_sample <- rast_sample[sel,]
  }
  rast_sample <- clean_NAsInf(rast_sample)
  # update attribute table to eliminate plots which include no information
  AttributeTable <- AttributeTable[unique(rast_sample$ID),]
  # get plot size
  nbPix_per_plot <- data.frame(table(rast_sample$ID))
  # only get common plots between nbPix_per_plot and nbPix_per_plot_init
  if (length(nbPix_per_plot$Freq)>0){
    whichPlots2keep <- as.numeric(intersect(nbPix_per_plot[,'Var1'],
                                            nbPix_per_plot_init[,'Var1']))

    nbPix_per_plot <- nbPix_per_plot[which(nbPix_per_plot[,'Var1']%in%whichPlots2keep),]
    nbPix_per_plot_init <- nbPix_per_plot_init[which(nbPix_per_plot_init[,'Var1']%in%whichPlots2keep),]
    Sunlit_plots <- which(nbPix_per_plot$Freq/nbPix_per_plot_init$Freq >= MinSun)
    nbPix_per_plot <- nbPix_per_plot[Sunlit_plots,]
    nbPix_per_plot_init <- nbPix_per_plot_init[Sunlit_plots,]
    selpix <- which(rast_sample$ID%in%nbPix_per_plot_init[,'Var1'])
    rast_sample <- rast_sample[selpix,]
  }
  return(list('rast_sample_vect' = rast_sample,
              'AttributeTable' = AttributeTable))
}
