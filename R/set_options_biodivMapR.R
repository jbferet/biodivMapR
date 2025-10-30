set_options <- function(fun, options = NULL){

  if (fun == 'biodivMapR_full_tiles'){
    if (is.null(options$nb_clusters))
      options$nb_clusters <- 50
    if (is.null(options$nb_samples_alpha))
      options$nb_samples_alpha <- 1e5
    if (is.null(options$nb_samples_beta))
      options$nb_samples_beta <- 2e3
    if (is.null(options$alphametrics))
      options$alphametrics <- 'shannon'
    if (is.null(options$Hill_order))
      options$Hill_order <- 1
    if (is.null(options$FDmetric))
      options$FDmetric <- NULL
    if (is.null(options$nb_iter))
      options$nb_iter <- 10
    if (is.null(options$nb_iter))
      options$pcelim <- 0.02
    if (is.null(options$maxRows))
      options$maxRows <- 1200
    if (is.null(options$moving_window))
      options$moving_window <- FALSE
    if (is.null(options$mosaic_output))
      options$mosaic_output <- TRUE
    if (is.null(options$weightIRQ))
      options$weightIRQ <- 4
  }
  return(options)
}




