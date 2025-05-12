#' extract pixel information from a raster based on SpatVectorCollection
#'
#' @param SpatVector SpatVectorCollection object
#' @param input_rast SpatRaster input raster
#' @param input_mask SpatRaster optional mask
#' @param min_sun numeric. minimum sunlit/unmasked proportion
#' @param prog boolean progressbar ?
#'
#' @return list including rast_sample: dataframe corresponding to
#' information extracted from input_rast, and AttributeTable from SpatVector

#' @importFrom terra sources extract
#' @export

extract_svc_from_rast <- function(SpatVector, input_rast,
                                  input_mask = NULL, min_sun = 0.25,
                                  prog = TRUE){

  # extract from list of SpatVectors in collection
  rastext <- lapply(X = SpatVector,
                    FUN = extract_vect_from_rast,
                    input_rast = input_rast,
                    input_mask = input_mask,
                    min_sun = min_sun, prog = prog)
  # update plot ID in collection
  nbPlots_total <- 0
  for (ind_vect in seq_len(length(SpatVector))){
    AttributeTable <- rastext[[ind_vect]]$AttributeTable
    rast_sample <- rastext[[ind_vect]]$rast_sample
    AttributeTable$ID_biodivMapR <- AttributeTable$ID_biodivMapR + nbPlots_total
    rast_sample$ID <- rast_sample$ID + nbPlots_total
    nbPlots_total <- max(AttributeTable$ID_biodivMapR)
    rastext[[ind_vect]]$AttributeTable <- AttributeTable
    rastext[[ind_vect]]$rast_sample <- rast_sample
  }
  rast_sample <- lapply(rastext,'[[','rast_sample')
  AttributeTable <- lapply(rastext,'[[','AttributeTable')
  rast_sample <- do.call(rbind,rast_sample)
  Attributes <- do.call(rbind,AttributeTable)
  return(list('rast_sample_vect' = rast_sample,
              'AttributeTable' = AttributeTable))
}
