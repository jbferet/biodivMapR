#' computes diversity metrics from a classification raster containing information
#' on a specific class for each pixel
#'
#' @param input_raster_path character. path for the input rasters
#' @param output_dir character. path for the output files
#' @param window_size numeric. window size for square plots
#' @param input_mask_path character. path for mask file
#' @param site_name character. nname of site
#' @param alpha_metrics list. alpha diversity metrics
#' @param beta_metrics list. beta diversity metrics
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param compute_beta boolean. set TRUE if beta is to be computed
#' @param nb_samples_beta numeric. number of samples to compute beta diversity
#'
#' @return diversity_maps_ground
#' @importFrom terra vect crs rast extract values aggregate writeRaster
#' @export

biodivMapR_classif <- function(
    input_raster_path, output_dir, window_size, input_mask_path = NULL, 
    site_name = NULL, alpha_metrics = c('richness', 'shannon', 'simpson', 'hill'), 
    beta_metrics = c('bray', 'brayturn', 'simpson_diss', 'sorensen', 'jaccard', 'jaccardturn'),
    pcelim = 0.02, compute_beta = TRUE, nb_samples_beta = 1000){

  # create output_dir if does not exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # read classification map
  input_rast <- terra::rast(input_raster_path)
  nb_clusters <- length(na.omit((unique(c(terra::values(input_rast))))))

  # produce a plot grid over the full raster
  grid_plot <- define_grid(raster_path = input_raster_path,
                           cellsize = window_size*terra::res(input_rast)[1])
  plots <- grid_plot$plots

  # account for mask
  input_mask <- NULL
  if (!is.null(input_mask_path))
    input_mask <- terra::rast(input_mask_path)

  # define plot selection for beta diversity and sample from raster
  Beta_info <- NULL
  if (compute_beta){
    if (nb_samples_beta > length(plots))
      nb_samples_beta <- length(plots)

    plots_beta <- plots[sample(x = seq_along(plots), nb_samples_beta,
                               replace = FALSE)]
    samples <- lapply(X = plots_beta, FUN = get_samples_from_plots,
                      y = terra::rast(input_raster_path))
    
    ll <- unlist(lapply(samples, length))
    sel <- which(ll>0.5*window_size*window_size)
    samples <- samples[sel]

    # compute spectral dissimilarity
    ssd <- lapply(X = samples,FUN = table)
    ssd <- lapply(X = ssd,FUN = get_normalized_ssd,
                  nb_clusters = nb_clusters, pcelim = pcelim)
    ssd <- do.call(rbind,ssd)
    mat_diss <- compute_dissimilarity(A = ssd, B = ssd, beta_metrics = beta_metrics)
    beta_pco <- list()
    for (beta in beta_metrics){
      mat_diss_dist <- stats::as.dist(mat_diss[[beta]], diag = FALSE, upper = FALSE)
      beta_pco[[beta]] <- pco(mat_diss_dist, k = 3)
    }
    Beta_info <- list('ssd' = ssd, 'mat_diss' = mat_diss, 'beta_pco' = beta_pco)
    # # save spectral dissimilarity
    # if (is.null(Beta_info_save))
    #   Beta_info_save <- file.path(output_dir, 'Beta_info_classif.RData')
    # save(Beta_info, file = Beta_info_save)
  }

  # compute alpha and beta diversity over the full image
  input_rast_tmp <- input_rast
  terra::values(input_rast_tmp) <- 0
  output_rast_tmp <- terra::aggregate(x = input_rast_tmp, fact = window_size)
  # get the id of the cells corresponding to the plots
  get_cell_plot <- function(x, y){
    x <- terra::vect(x)
    x <- terra::project(x = x, y = y)
    res <- terra::extract(x = y, y = x, cell = TRUE)
    res <- c(unlist(res[1,]))
    return(res)
  }
  extracted_val <- lapply(X = plots, FUN = get_cell_plot,
                          y = output_rast_tmp)

  # get alpha and beta diversity metrics
  samples <- lapply(X = plots, FUN = get_samples_from_plots,
                    y = terra::rast(input_raster_path))

  alphabeta <- alphabeta_window_classif(SSwindow = samples,
                                        nb_clusters = nb_clusters,
                                        Beta_info = Beta_info,
                                        alpha_metrics = alpha_metrics,
                                        beta_metrics = beta_metrics,
                                        Hill_order = 1, pcelim = pcelim)
  cell_order <- unlist(lapply(extracted_val, '[[', 'cell'))
  if (is.null(Beta_info)){
    elim <- grep(x = names(alphabeta), pattern = 'pcoa_')
    if (length(elim) > 0)
      alphabeta <- alphabeta[-elim]
  }


  # save diversity maps
  diversity_maps_ground <- list()
  for (idx in names(alphabeta)){
    if (length(grep(x = idx, pattern = 'pcoa_'))>0 & !is.null(Beta_info)){
      beta <- list(output_rast_tmp, output_rast_tmp, output_rast_tmp)
      # for (i in 1:3)
      #   beta[[i]][cell_order] <- unlist(lapply(alphabeta[[idx]], '[[', i))
      for (i in 1:3)
        beta[[i]][cell_order] <- alphabeta[[idx]][,i]
      beta <- terra::rast(beta)
      names(beta) <- c('pco1', 'pco2', 'pco3')
      filename <- file.path(output_dir, paste0(idx,'_classif.tiff'))
      if (!is.null(site_name))
        filename <- file.path(output_dir,
                              paste0(site_name, '_',idx,'_classif.tiff'))
      terra::writeRaster(x = beta, filename = filename, overwrite = TRUE)
      diversity_maps_ground[[idx]] <- filename
    } else {
      alpha <- output_rast_tmp
      alpha[cell_order] <-alphabeta[[idx]]
      filename <- file.path(output_dir, paste0(idx, '_classif.tiff'))
      names(alpha) <- idx
      if (!is.null(site_name))
        filename <- file.path(output_dir, paste(site_name, idx,
                                                'classif.tiff', sep = '_'))
      terra::writeRaster(x = alpha, filename = filename, overwrite = TRUE)
      diversity_maps_ground[[idx]] <- filename
    }
  }
  return(diversity_maps_ground)
}
