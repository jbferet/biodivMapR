#' computes diversity metrics from raster
#'
#' @param input_raster_path character. path for the input rasters
#' @param output_dir character. path for the output files
#' @param window_size numeric. window size for square plots
#' @param input_mask_path character. path for mask file
#' @param site_name character. nname of site
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param compute_beta boolean. set TRUE if beta is to be computed
#' @param nb_samples_beta numeric. number of samples to compute beta diversity
#'
#' @return diversity_maps_ground
#' @importFrom dissUtils diss
#' @importFrom terra vect crs rast extract values aggregate writeRaster
#' @export

biodivMapR_full_classif <- function(input_raster_path, output_dir, window_size,
                                    input_mask_path = NULL, site_name = NULL,
                                    pcelim = 0.02, compute_beta = TRUE,
                                    nb_samples_beta = 1000){

  # define all alpha metrics
  alpha_metrics <- c('richness', 'shannon', 'simpson', 'hill')
  input_rast <- terra::rast(input_raster_path)
  nb_clusters <- length((unique(c(terra::values(input_rast)))))

  # produce a plot grid over the full raster
  grid_plot <- define_grid(raster_path = input_raster_path,
                           cellsize = window_size)
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

    # compute spectral dissimilarity
    ssd <- lapply(X = samples,FUN = table)
    ssd <- lapply(X = ssd,FUN = get_normalized_ssd,
                  nb_clusters = nb_clusters, pcelim = pcelim)
    ssd <- do.call(rbind,ssd)
    mat_bc <- dissUtils::diss(ssd, ssd, method = 'braycurtis')
    # ssd_list <- list(ssd, ssd)
    # mat_bc <- compute_bc_diss(ssd_list, pcelim)
    Beta_info <- list('SSD' = ssd, 'MatBC' = mat_bc)
    mat_bc_dist <- stats::as.dist(mat_bc, diag = FALSE, upper = FALSE)
    BetaPCO <- pco(mat_bc_dist, k = 3)
    Beta_info <- list('SSD' = ssd,
                      'MatBC' = mat_bc,
                      'BetaPCO' = BetaPCO)
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
    terra::crs(x) <- terra::crs(y)
    res <- terra::extract(x = y, y = x, cell = TRUE)
    res <- c(unlist(res))
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
                                        Hill_order = 1, pcelim = pcelim)
  cell_order <- unlist(lapply(extracted_val, '[[', 'cell'))
  if (is.null(Beta_info)){
    elim <- which(names(alphabeta)=='PCoA_BC')
    if (length(elim) > 0)
      alphabeta <- alphabeta[-elim]
  }


  # save diversity maps
  diversity_maps_ground <- list()
  for (idx in names(alphabeta)){
    if (idx == 'PCoA_BC' & !is.null(Beta_info)){
      beta <- list(output_rast_tmp, output_rast_tmp, output_rast_tmp)
      for (i in 1:3)
        beta[[i]][cell_order] <- unlist(lapply(alphabeta[[idx]], '[[', i))
      beta <- terra::rast(beta)
      names(beta) <- c('pco1', 'pco2', 'pco3')
      filename <- file.path(output_dir, 'beta_classif.tiff')
      if (!is.null(site_name))
        filename <- file.path(output_dir,
                              paste0(site_name, '_beta_classif.tiff'))
      terra::writeRaster(x = beta, filename = filename, overwrite = TRUE)
      diversity_maps_ground$beta <- filename
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
