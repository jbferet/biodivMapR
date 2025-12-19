#' computes diversity metrics from raster
#'
#' @param input_raster_path character. path for the input rasters
#' @param output_dir character. path for the output files
#' @param input_mask_path character. path for mask file
#' @param site_name character. nname of site
#' @param alpha_metrics character.
#' @param nb_samples_beta numeric. number of samples to compute beta diversity
#'
#' @return diversity_maps_ground
#' @importFrom dissUtils diss
#' @export

biodivMapR_sdm <- function(input_raster_path, output_dir,
                           input_mask_path = NULL, site_name = NULL,
                           alpha_metrics = 'shannon',
                           nb_samples_beta = 1000){

  # define all alpha metrics
  input_rast <- terra::rast(input_raster_path)
  nb_clusters <- dim(input_rast)[3]

  # account for mask
  input_mask <- 1+0*input_rast[[1]]
  if (!is.null(input_mask_path))
    input_mask <- terra::rast(input_mask_path)

  # define plot selection for beta diversity and sample from raster
  points <-terra::spatSample(input_mask, nb_samples_beta, as.points=TRUE)
  # check how many data points and adjust sampling
  tab <- table(terra::values(points))
  ratio <- tab['1']/nb_samples_beta
  if (ratio<0.95){
    nb_samples_beta_adjust <- 1.2*nb_samples_beta/ratio
    points <-terra::spatSample(input_mask, nb_samples_beta_adjust, as.points = TRUE)
    points <- points[which(terra::values(points) ==1)]
    # check how many data points and adjust sampling
    tab <- table(terra::values(points))
    if (tab['1']>nb_samples_beta)
      points <- points[1:nb_samples_beta]
  }
  ssd <- terra::extract(x = input_rast, y = points)
  ssd$ID <- NULL
  names(ssd) <- seq_len(nb_clusters)
  ssd <- as.matrix(ssd)
  # compute spectral dissimilarity
  mat_bc <- dissUtils::diss(ssd, ssd, method = 'braycurtis')
  # ssd_list <- list(ssd, ssd)
  # mat_bc <- compute_bc_diss(ssd_list, pcelim = 0)
  mat_bc_dist <- stats::as.dist(mat_bc, diag = FALSE, upper = FALSE)
  BetaPCO <- pco(mat_bc_dist, k = 3)
  Beta_info <- list('SSD' = ssd, 'MatBC' = mat_bc, 'BetaPCO' = BetaPCO)

  # compute alpha and beta diversity over the full image
  # terra::values(output_rast_tmp) <- 0
  extracted_val <- terra::values(input_rast)
  extracted_mask <- terra::values(input_mask)
  extracted_val <- extracted_val[which(extracted_mask==1), ]
  colnames(extracted_val) <- seq_len(nb_clusters)
  extracted_val <- lapply(seq_len(nrow(extracted_val)),
                          function(i) extracted_val[i,])

  alphabeta <- alphabeta_window_sdm(ssd = extracted_val,
                                    Beta_info = Beta_info,
                                    alpha_metrics = alpha_metrics,
                                    Hill_order = 1)
  cell_order <- which(extracted_mask==1)

  # save diversity maps
  diversity_maps_ground <- list()
  output_rast_tmp <- NA*input_rast[[1]]
  for (idx in names(alphabeta)){
    if (idx == 'PCoA_BC'){
      beta <- list(output_rast_tmp, output_rast_tmp, output_rast_tmp)
      for (i in 1:3)
        beta[[i]][cell_order] <- unlist(lapply(alphabeta[[idx]], '[[', i))
      beta <- terra::rast(beta)
      names(beta) <- c('pco1', 'pco2', 'pco3')
      filename <- file.path(output_dir, 'beta_sdm.tiff')
      if (!is.null(site_name))
        filename <- file.path(output_dir, paste0(site_name, '_beta_sdm.tiff'))
      terra::writeRaster(x = beta, filename = filename, overwrite = TRUE)
      diversity_maps_ground$beta <- filename
    } else {
      alpha <- output_rast_tmp
      alpha[cell_order] <-alphabeta[[idx]]
      filename <- file.path(output_dir, paste0(idx, '_sdm.tiff'))
      names(alpha) <- idx
      if (!is.null(site_name))
        filename <- file.path(output_dir, paste(site_name, idx,
                                                'sdm.tiff', sep = '_'))
      terra::writeRaster(x = alpha, filename = filename, overwrite = TRUE)
      diversity_maps_ground[[idx]] <- filename
    }
  }
  return(diversity_maps_ground)
}
