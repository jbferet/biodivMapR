#' get spectral species corresponding to polygons in a SpatVector object
#
#' @param SpatVector SpatVector.
#' @param input_rast SpatRaster.
#' @param Kmeans_info list. kmeans description obtained from function get_kmeans
#' @param SelectBands numeric. bands selected from input data
#' @param Functional boolean. should functional diversity be computed as well?
#' @param input_mask SpatRaster.
#' @param rast_sample dataframe.
#' @param AttributeTable dataframe.
#' @param MinSun numeric. minimum amount of sunlit pixels in the plots
#
#' @return list
#' @importFrom dplyr group_by
#' @importFrom tidyr nest
#' @importFrom terra extract sources values
#' @export
#'
spectralspecies_per_polygon <- function(SpatVector, input_rast,
                                        Kmeans_info, SelectBands = NULL,
                                        Functional = NULL, input_mask = NULL,
                                        rast_sample = NULL, AttributeTable = NULL,
                                        MinSun = 0.25){

  FRic <- FEve <- FDiv <- FunctDiv <- NULL
  SSValid <- NULL
  # extract pixel info from vector data
  if (is.null(rast_sample)){
    rastext <- extract_vect_from_rast(SpatVector = SpatVector,
                                      input_rast = input_rast,
                                      input_mask = input_mask,
                                      MinSun = MinSun,
                                      prog = F)
    # update plot ID in collection
    rast_sample <- rastext$rast_sample_vect
    AttributeTable <- rastext$AttributeTable
  } else if (is.null(rast_sample) | is.null(AttributeTable)){
    warning('"rast_sample" or "AttributeTable" missing for "spectralspecies_per_polygon"')
  }

  nbPix_per_plot <- data.frame(table(rast_sample$ID))
  # only get common plots between nbPix_per_plot and nbPix_per_plot_init
  if (length(nbPix_per_plot$Freq)>0){
    if (is.null(SelectBands)) SelectBands <- seq_len(dim(rast_sample)[2]-1)
    rast_sample_noID <- rast_sample
    rast_sample_noID$ID <- NULL
    SSValid <- get_spectralSpecies(inputdata = rast_sample_noID,
                                   Kmeans_info = Kmeans_info,
                                   SelectBands = SelectBands)
    SSValid$win_ID <- rast_sample$ID
    # Functional diversity
    if (!is.null(Functional)){
      if (is.null(SelectBands)) SelectBands <- seq_len(ncol(rast_sample_noID))
      # center reduce data
      inputdata_cr <- center_reduce(X = rast_sample_noID[SelectBands],
                                    m = Kmeans_info$MinVal,
                                    sig = Kmeans_info$Range)
      inputdata_cr$ID <- rast_sample$ID
      inputdata_cr <- inputdata_cr %>% split(.$ID)
      inputdata_cr <- lapply(inputdata_cr, function(x) x[ , !(names(x) %in% "ID")])
      # in case only one dimension
      inputdata_cr <- lapply(inputdata_cr, data.frame)
      FunctDiv <- lapply(X = inputdata_cr,
                         FUN = get_functional_diversity,
                         FDmetric = Functional)
      # FunctDiv <- data.frame('FRic' = unlist(lapply(FunctDiv, '[[','FRic')),
      #                        'FEve' = unlist(lapply(FunctDiv, '[[','FEve')),
      #                        'FDiv' = unlist(lapply(FunctDiv, '[[','FDiv')))
      FunctDiv <- data.frame('FRic' = unlist(lapply(FunctDiv, '[[',1)),
                             'FEve' = unlist(lapply(FunctDiv, '[[',2)),
                             'FDiv' = unlist(lapply(FunctDiv, '[[',3)),
                             'FDis' = unlist(lapply(FunctDiv, '[[',4)),
                             'FRaoq' = unlist(lapply(FunctDiv, '[[',5)))
    }
  } else if (inherits(SpatVector, what = 'SpatVector')){
    AttributeTable <- values(SpatVector)
    AttributeTable$source <- basename(sources(SpatVector))
  }
  return(list('SSValid' = SSValid,
              'AttributeTable' = AttributeTable,
              'FunctDiv' = FunctDiv))
}
