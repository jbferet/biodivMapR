#' extract informaton from raster to prepare for SFS
#'
#' @param input_raster SpatRaster or list of SpatRaster
#' @param input_mask SpatRaster corresponding to mask
#' @param min_sun numeric.
#' @param obs_vect SpatVector or SpatVectorCollection
#'
#' @return list
#'
#' @export

extract_from_rast_sfs <- function(input_raster, input_mask, min_sun, obs_vect){

  # 2- extract information from field plots
  if (inherits(obs_vect, what = 'SpatVectorCollection')){
    # extract from list of SpatVectors in collection
    rastext <- lapply(X = obs_vect,
                      FUN = extract_vect_from_rast,
                      input_rast = input_raster,
                      input_mask = input_mask,
                      min_sun = min_sun, prog = FALSE)
    # update plot ID in collection
    nbPlots_total <- 0
    for (ind_vect in seq_len(length(obs_vect))){
      AttributeTable <- rastext[[ind_vect]]$AttributeTable
      rast_sample_vect <- rastext[[ind_vect]]$rast_sample_vect
      AttributeTable$ID_biodivMapR <- AttributeTable$ID_biodivMapR +
        nbPlots_total
      rast_sample_vect$ID <- rast_sample_vect$ID + nbPlots_total
      nbPlots_total <- max(AttributeTable$ID_biodivMapR)
      rastext[[ind_vect]]$AttributeTable <- AttributeTable
      rastext[[ind_vect]]$rast_sample_vect <- rast_sample_vect
    }
    rast_sample_vect <- lapply(rastext,'[[','rast_sample_vect')
    AttributeTable <- lapply(rastext,'[[','AttributeTable')
    rast_val <- do.call(rbind,rast_sample_vect)
    Attributes <- do.call(rbind,AttributeTable)
  } else if (inherits(obs_vect, what = 'SpatVector')){
    # extract from list of SpatVectors in collection
    rastext <- extract_vect_from_rast(SpatVector = obs_vect,
                                      input_rast = input_raster,
                                      input_mask = input_mask,
                                      min_sun = min_sun)
    # update plot ID in collection
    # rast_sample_vect <- lapply(rastext,'[[','rast_sample_vect')
    # AttributeTable <- lapply(rastext,'[[','AttributeTable')
    # rast_sample_vect <- rastext$rast_sample_vect
    # AttributeTable <- rastext$AttributeTable
    rast_val <- rastext$rast_sample_vect
    Attributes <- rastext$AttributeTable
    nbPlots_total <- length(Attributes$ID_biodivMapR)
  }
  IDplot <- rast_val$ID
  rast_val$ID <- NULL
  return(list('rast_val' = rast_val,
              'Attributes' = Attributes,
              'nbPlots_total' = nbPlots_total,
              'IDplot' = IDplot))
}
