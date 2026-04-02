#' explore performances for a range of cluster number
#'
#' @param repet character identifier of the repetition
#' @param input_raster_path SpatRaster or list of SpatRaster
#' @param input_vector_path SpatVector or SpatVectorCollection
#' @param input_mask_path SpatRaster corresponding to mask
#' @param selected_bands numeric. bnds to select from input_raster
#' @param obs_criterion character. richness, shannon, simpson or BC
#' @param alpha_metrics character. richness, shannon, simpson or BC
#' @param outputdir character. output directory
#' @param nbClust_list numeric.
#' @param min_sun numeric.
#' @param nb_iter numeric.
#' @param pcelim numeric.
#' @param nb_samples_alpha numeric.
#' @param Hill_order numeric.
#' @param algorithm character.
#' @param p list
#'
#' @return list divIndex_est

#' @importFrom terra vect rast
#' @export

explore_cluster_range <- function(repet, input_raster_path, input_vector_path,
                                  input_mask_path, selected_bands,
                                  obs_criterion = 'shannon', alpha_metrics,
                                  outputdir = './', nbClust_list, min_sun = 0.25,
                                  nb_iter = 10, pcelim = 0.02, nb_samples_alpha = 1e5,
                                  Hill_order = 1, algorithm = 'Hartigan-Wong',
                                  p = NULL){

  input_vector <- terra::vect(input_vector_path)
  input_rast <- terra::rast(x = input_raster_path)
  input_mask <- NULL
  if (!is.null(input_mask_path))
    input_mask <- terra::rast(x = input_mask_path)

  # extract information from SpatVectorCollection or SpatVector
  if (inherits(input_vector, what = 'SpatVectorCollection')){
    rastext <- extract_svc_from_rast(SpatVector = input_vector,
                                     input_rast = input_rast,
                                     input_mask = input_mask,
                                     min_sun = min_sun, prog = FALSE)
  } else if (inherits(input_vector, what = 'SpatVector')){
    rastext <- extract_vect_from_rast(SpatVector = input_vector,
                                      input_rast = input_rast,
                                      input_mask = input_mask,
                                      min_sun = min_sun, prog = FALSE)
  }
  rast_sample <- rastext$rast_sample_vect
  AttributeTable <- rastext$AttributeTable

  # perform kmeans for the full list of nbClust_list
  Kmeans_info <- explore_kmeans(input_rast = input_rast,
                                input_mask = input_mask,
                                selected_bands = selected_bands,
                                nbClust_list = nbClust_list,
                                nb_iter = nb_iter,
                                nb_samples_alpha = nb_samples_alpha,
                                algorithm = algorithm,
                                nbCPU = 1)
  Beta_info <- NULL
  divPlots_kmeans <- lapply(X = Kmeans_info,
                            FUN = get_diversity_from_plots_cluster_v2,
                            rast_sample = rast_sample,
                            AttributeTable = AttributeTable,
                            Hill_order = Hill_order,
                            Beta_info = Beta_info,
                            selected_bands = selected_bands,
                            alpha_metrics = alpha_metrics,
                            pcelim = pcelim)
  gc()

  divIndex_est <- list()
  for (crit0 in obs_criterion){
    Sel1 <- lapply(X = divPlots_kmeans, FUN = '[[', 'specdiv')
    if (crit0 %in% c('richness', 'shannon', 'simpson', 'hill')){
      divIndex_est[[crit0]] <- data.frame(lapply(X = Sel1,
                                                 FUN = '[[',
                                                 paste0(crit0, '_mean')))
      colnames(divIndex_est[[crit0]]) <- unlist(nbClust_list)
      # } else if (crit0 %in% 'BC'){
      #   divIndex_est[[crit0]] <- lapply(X = divPlots_kmeans,
      #                                            FUN = '[[', 'BC_dissimilarity')
      #   divIndex_est[[crit0]] <- lapply(X = divIndex_est[[crit0]],
      #                                            FUN = 'as.dist')
      #   divIndex_est[[crit0]] <- data.frame(lapply(X = divIndex_est[[crit0]],
      #                                                       FUN = 'c'))
      #   colnames(divIndex_est[[crit0]]) <- unlist(nbClust_list)
    }
    # save results to ease analysis if stopped at some point
    fileName <- file.path(outputdir, paste0(crit0,'_',repet,'.csv'))
    readr::write_delim(x = round(divIndex_est[[crit0]], digits = 5),
                       file = fileName, delim = '\t', progress = FALSE)
  }
  if (!is.null(p))
    p()
  rm(list = ls(all=TRUE)[-which(ls()=='divIndex_est')])
  gc()
  return(divIndex_est)

}
