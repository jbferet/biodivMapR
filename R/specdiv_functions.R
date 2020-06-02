# ==============================================================================
# biodivMapR
# specdiv_functions.R
# ==============================================================================
#### Spectral diversity functions
# Etienne Laliberte, February 2019
# original script can be found here: https://github.com/elaliberte/specdiv
# ==============================================================================
# This Library performs partitioning of plant spectral diversity into alpha and
# beta components
# it follows the method proposed by Lalibert?, Schweiger & Legendre (2020),
# Partitioning plant spectral diversity into alpha and beta components,
# Ecology letters (https://doi.org/10.1111/ele.13429)
# ==============================================================================



#' Brightness normalization -----
#'
#' @param x numeric. spectrum to be normalized
#'
#' @return normalized spectrum
#' @export

bright_norm <- function(x) {
  x_norm <- x / sqrt(base::sum(x^2))
  return(x_norm)
}


#' Brightness normalization for data frame ----
#'
#' @param x dataframe. spectrum to be normalized
#'
#' @return dataframe including normalized spectrum
#' @export

bright_norm_df <- function(x) {
  x$refl_norm <- x$refl / sqrt(base::sum(x$refl^2))
  return(x)
}


#' #' Short PCA for matrix-----
#' #'
#' #' @param x numeric. matrix on which PCA should be applied
#' #' @param scaling numeric.
#' #' @param p numeric. cumulated proportion of variance explained by the components provided as outputs
#' #'
#' #' @return out list. list of PCA outputs
#' #' @export
#'
#' pca_mat <- function(x, scaling = c(1, 2), p = 0.99) {
#'   Y <- scale(x, center = TRUE, scale = FALSE)
#'   n <- nrow(Y)
#'   Y.svd = svd(Y)
#'   values = (1/(n - 1)) * Y.svd$d^2
#'   epsilon = sqrt(.Machine$double.eps)
#'   k <- base::sum(values > epsilon)
#'   values <- values[1:k]
#'   prop <- values / base::sum(values)
#'   cumprop = cumsum(prop)
#'   if (p < cumprop[1]) which.values <- c(1, 2) else which.values <- which(cumprop < p)
#'   values.sel <- values[c(which.values, length(which.values) + 1)]
#'   n.pcs <- length(values.sel)
#'   U <- as.matrix(Y.svd$v[, 1:n.pcs])
#'   if (scaling == 1) {
#'     obj <- Y %*% U
#'     descript <- U
#'   }
#'   else {
#'     obj <- sqrt(n - 1) * as.matrix(Y.svd$u[, 1:n.pcs])
#'     descript <- U %*% diag(values.sel^(0.5))
#'   }
#'   colnames(obj) <- paste0('PC', 1:n.pcs)
#'   colnames(descript) <- paste0('PC', 1:n.pcs)
#'   rownames(descript) <- colnames(x)
#'   prop.sel <- prop[1:n.pcs]; names(prop.sel) <- colnames(descript)
#'   cumprop.sel <- cumprop[1:n.pcs]; names(cumprop.sel) <- colnames(descript)
#'   out <- list(obj = obj, descript = descript, prop = prop.sel, cumprop = cumprop.sel)
#'   return(out)
#' }



#' #' PCA for image cube -----
#' #'
#' #' @param cube numeric. matrix on which PCA should be applied
#' #' @param scaling numeric.
#' #' @param p numeric. cumulated proportion of variance explained by the components provided as outputs
#' #'
#' #' @return out list. list of PCA outputs
#' #' @export
#' pca <- function(cube, scaling = c(1, 2), p = 0.99) {
#'   require(tidyverse)
#'   require(raster)
#'   require(sp)
#'   points <- rasterToPoints(cube, spatial = F)
#'   xy <- points[, 1:2]
#'   pixels <- points[, 3:ncol(points)]
#'   Y <- scale(pixels, center = TRUE, scale = FALSE)
#'   n <- nrow(Y)
#'   Y.svd = svd(Y)
#'   values = (1/(n - 1)) * Y.svd$d^2
#'   epsilon = sqrt(.Machine$double.eps)
#'   k <- base::sum(values > epsilon)
#'   values <- values[1:k]
#'   prop <- values / base::sum(values)
#'   cumprop = cumsum(prop)
#'   if (p < cumprop[1]) which.values <- c(1, 2) else which.values <- which(cumprop < p)
#'   values.sel <- values[c(which.values, length(which.values) + 1)]
#'   n.pcs <- length(values.sel)
#'   U <- as.matrix(Y.svd$v[, 1:n.pcs])
#'   if (scaling == 1) {
#'     obj <- Y %*% U
#'     descript <- U
#'   }
#'   else {
#'     obj <- sqrt(n - 1) * as.matrix(Y.svd$u[, 1:n.pcs])
#'     descript <- U %*% diag(values.sel^(0.5))
#'   }
#'   colnames(obj) <- paste0('PC', 1:n.pcs)
#'   colnames(descript) <- paste0('PC', 1:n.pcs)
#'   rownames(descript) <- colnames(pixels)
#'   points_df <- SpatialPixelsDataFrame(xy, as.data.frame(obj), proj4string = CRS(proj4string(cube) ))
#'   prop.sel <- prop[1:n.pcs]; names(prop.sel) <- colnames(descript)
#'   cumprop.sel <- cumprop[1:n.pcs]; names(cumprop.sel) <- colnames(descript)
#'   cube_pc <- brick(points_df)
#'   out <- list(cube_pc = cube_pc, band_contrib = descript, prop = prop.sel, cumprop = cumprop.sel)
#'   return(out)
#' }


# ### Make PC plots ----
# make_plot <- function(df) {
#   x <- ggplot(df, aes(x = x, y = y) ) +
#     geom_raster(aes(fill = value)) +
#     scale_fill_gradientn(colors = rainbow(20)) +
#     ggtitle(label = unique(df$PC)) +
#     theme_void() +
#     theme(legend.position = 'none',
#           plot.title = element_text(face = 'bold', hjust = 0.5, size = 15) ) +
#     coord_equal()
#   return(x)
# }


# ### Build the PC plots ----
# build_pc_plot <- function(cube) {
#   require(tidyverse)
#   points <- rasterToPoints(cube, spatial = F) %>%
#     as_tibble() %>%
#     gather(key = PC, value = value, -x, -y) %>%
#     mutate(PC = factor(PC, levels = paste0('PC', 1:n_distinct(PC)) ) )
#   cube_plots <- points %>%
#     group_by(PC) %>%
#     do(plots = make_plot(df = .))
# return(cube_plots)
# }


# ### Show the PC plots ----
# show_pc_plot <- function(x, ...) {
#   require(ggplot2)
#   require(gridExtra)
#   grid.arrange(grobs = x$plots, ...)
# }



#' Count non-NA rows
#'
#' @param x dataframe. includes a list of values
#'
#' @return dataframe number of non-NA values
#' @importFrom stats na.omit
#' @export

count_noNA <- function(x) {
  n <- nrow(na.omit(x))
  n <- as.data.frame(n)
}


#' Count masked/missing pixels -----
#'
#' @param cube dataframe. data cube
#' @param fact numeric. window size
#'
#' @return cube_count
#' @import raster
#' @import tidyverse
#' @export

count_pixels <- function(cube, fact = 40) {
  require(raster)
  require(tidyverse)
  # Get plots (or communities)
  new_res <- res(cube) * fact
  cube_plots <- raster(crs = proj4string(cube))
  extent(cube_plots) <- extent(cube)
  res(cube_plots) <- new_res
  cube_plots <- setValues(cube_plots, values = 1:ncell(cube_plots))
  plot_xy <- rasterToPoints(cube_plots)[,1:2]
  cube_pixels <- disaggregate(cube_plots, fact = fact)
  # Convert to points
  plot_points <- rasterToPoints(cube_pixels) %>%
    as_tibble() %>%
    rename(group = layer)
  cube_points <- rasterToPoints(cube) %>%
    as_tibble() %>%
   right_join(plot_points, by = c('x', 'y'))
  n_pixels <- cube_points %>%
    group_by(group) %>%
    do(count_noNA(.)) %>%
    mutate(n_total = fact^2,
           prop = n / n_total) %>%
    ungroup()
  group_df <- SpatialPixelsDataFrame(plot_xy, dplyr::select(n_pixels, n:prop), proj4string = CRS(proj4string(cube) ))
  cube_count <- brick(group_df)
  return(cube_count)
}

#' SS gamma and alpha ----
#'
#' @param Y numeric. data matrix
#'
#' @return out
#' @export

sum_squares <- function(Y) {
  n <- nrow(Y)
  Y.cent <- scale(Y, center = T, scale = F)
  sij <- Y.cent^2
  SS.total <- base::sum(sij)
  SS.row <- rowSums(sij)
  SS.col <- colSums(sij)
  fcsd <- SS.col / SS.total
  lcsd <- SS.row / SS.total
  sdiv <- SS.total / (n - 1)
  out <- list(ss = SS.total, sdiv = sdiv, lcsd = lcsd, fcsd = fcsd)
  return(out)
}


#' SS beta ----
#'
#' @param Y numeric. data matrix
#' @param m numeric. unused
#'
#' @return out
#' @export

sum_squares_beta <- function(Y, m) {
  n <- nrow(Y)
  Y.cent <- bind_cols(dplyr::select(Y, group), as.data.frame(scale(dplyr::select(Y, -group), scale = F)))
  mskj <- Y.cent %>%
    group_by(group) %>%
    mutate_at(vars(-group_cols()), function(x) (mean(x))^2) %>%
    summarise_at(vars(-group_cols()), sum) %>%
    ungroup() %>%
    dplyr::select(-group)
  SSbk <- rowSums(mskj)
  SSbj <- colSums(mskj)
  SSb <- base::sum(SSbk)
  sdiv <- SSb / (n - 1)
  fcsd <- SSbj / SSb
  lcsd <- SSbk / SSb
  out <- list(ss = SSb, sdiv = sdiv, lcss = SSbk, lcsd = lcsd, fcsd = fcsd)
  return(out)
}

#' Partitioning spectral diversity ----
#'
#' @param cube dataframe. data cube corresponding to raster
#' @param fact numeric. window size
#' @param prop numeric. proportion of valid pixels in a window under which window is discarded
#' @param n numeric.
#'
#' @return out list. includes rasters and mean alpha, beta and gamma
#' @import raster
#' @import tidyverse
#' @export

specdiv <- function(cube, fact = 40, prop = 0.5, n = 1) {
  require(raster)
  require(tidyverse)

  # Find minimum number for resampling
  n_pixels <- count_pixels(cube, fact = fact)

  # Remove plots with n pixels less than prop
  plot_mask <- subset(n_pixels, 'prop') < prop
  n_pixels_mask <- raster::mask(n_pixels, plot_mask, maskvalue = 1)
  min_pixels <- minValue(n_pixels_mask, 1)

  # Get cube with plots
  nlayers <- dim(cube)[3]
  cube_plots <- raster(crs = proj4string(n_pixels_mask))
  extent(cube_plots) <- extent(n_pixels_mask)
  res(cube_plots) <- res(n_pixels_mask)
  cube_plots <- setValues(cube_plots, values = 1:ncell(cube_plots))
  cube_plots_masked <- raster::mask(cube_plots, plot_mask, maskvalue = 1)
  cube_pixels <- disaggregate(cube_plots_masked, fact = fact)

  # Convert to points
  plots_points <- rasterToPoints(cube_plots_masked) %>%
    as_tibble() %>%
    rename(group = layer)
  pixels_points <- rasterToPoints(cube_pixels) %>%
    as_tibble() %>%
    rename(group = layer)

  # Objects to store results
  gamma_ss <- double()
  gamma_sdiv <- double()
  gamma_fcsd <- matrix(nrow = n, ncol = nlayers, dimnames = list(1:n, names(cube)))
  alpha_sdiv <- tibble()
  alpha_fcsd <- tibble()
  alpha_ss <- tibble()
  beta_ss <- double()
  beta_sdiv <- double()
  beta_fcsd <- matrix(nrow = n, ncol = nlayers, dimnames = list(1:n, names(cube)))
  beta_lcsd <- tibble()
  beta_lcss <- tibble()

  # Loop to randomly sample min pixels
  for (i in 1:n) {
    cube_points <- rasterToPoints(cube) %>%
      as_tibble() %>%
      inner_join(pixels_points, by = c('x', 'y')) %>%
      group_by(group) %>%
      sample_n(size = min_pixels) %>%
      ungroup()
    xy_all <- dplyr::select(cube_points, x, y)
    xy_plots <- dplyr::select(plots_points, x, y)

    # Gamma diversity
    cube_points_sel_gamma <- cube_points %>%
      dplyr::select(-x, -y, -group)
    sdiv_gamma <- sum_squares(cube_points_sel_gamma)
    gamma_ss[i] <- sdiv_gamma$ss
    gamma_sdiv[i] <- sdiv_gamma$sdiv
    gamma_fcsd[i, ] <- sdiv_gamma$fcsd

    # Alpha diversity
    cube_points_sel_alpha <- cube_points %>%
      dplyr::select(-x, -y)
    sdiv_alpha <- cube_points_sel_alpha %>%
      group_by(group) %>%
      do(res = sum_squares(Y = dplyr::select(., -group)) )
    # Get sdiv for each community
    alpha_sdiv_tmp <- tibble(rep = i, group = sdiv_alpha$group, sdiv = sapply(sdiv_alpha$res, function(x) x$sdiv))
    # store
    alpha_sdiv <- bind_rows(alpha_sdiv, alpha_sdiv_tmp)

    # Get fcsd for each community
    alpha_fcsd_tmp <- bind_cols(rep = rep(i, length(sdiv_alpha$group)), group = sdiv_alpha$group, as.data.frame(t(sapply(sdiv_alpha$res, function(x) x$fcsd))))
    alpha_fcsd <- bind_rows(alpha_fcsd, alpha_fcsd_tmp)

    # Get ss for each community
    alpha_ss_tmp <- tibble(rep = i, group = sdiv_alpha$group, ss = sapply(sdiv_alpha$res, function(x) x$ss))
    alpha_ss <- bind_rows(alpha_ss, alpha_ss_tmp)

    # Beta diversity
    cube_points_sel_beta <- cube_points %>%
      dplyr::select(-x, -y)
    sdiv_beta <- sum_squares_beta(cube_points_sel_beta, m = min_pixels)
    beta_ss[i] <- sdiv_beta$ss
    beta_sdiv[i] <- sdiv_beta$sdiv
    beta_fcsd[i, ] <- sdiv_beta$fcsd
    beta_lcsd_tmp <- tibble(rep = i, group = 1:length(sdiv_beta$lcsd), lcsd = sdiv_beta$lcsd)
    beta_lcsd <- bind_rows(beta_lcsd, beta_lcsd_tmp)
    beta_lcss_tmp <- tibble(rep = i, group = 1:length(sdiv_beta$lcss), lcss = sdiv_beta$lcss)
    beta_lcss <- bind_rows(beta_lcss, beta_lcss_tmp)
  }

  # Get results together

  # LCSD beta
  lcsd_beta_values <- beta_lcsd %>%
    rename(lcsd_beta = lcsd) %>%
    dplyr::select(-rep) %>%
    group_by(group) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    dplyr::select(-group)
  lcsd_beta_df <- SpatialPixelsDataFrame(xy_plots, lcsd_beta_values, proj4string = CRS(proj4string(cube) ))
  lcsd_beta_raster <- raster(lcsd_beta_df)

  # LCSS beta
  lcss_beta_values <- beta_lcss %>%
    rename(lcss_beta = lcss) %>%
    dplyr::select(-rep) %>%
    group_by(group) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    dplyr::select(-group)
  lcss_beta_df <- SpatialPixelsDataFrame(xy_plots, lcss_beta_values, proj4string = CRS(proj4string(cube) ))
  lcss_beta_raster <- raster(lcss_beta_df)

  # FCSD
  fcsd_beta <- colMeans(beta_fcsd)
  fcsd_gamma <- colMeans(gamma_fcsd)
  fcsd_alpha_values <- alpha_fcsd %>%
    dplyr::select(-rep) %>%
    group_by(group) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    dplyr::select(-group)
  fcsd_alpha_df <- SpatialPixelsDataFrame(xy_plots, fcsd_alpha_values, proj4string = CRS(proj4string(cube) ))
  fcsd_alpha_brick <- brick(fcsd_alpha_df)
  fcsd_alpha_mean <- colMeans(fcsd_alpha_values)

  # SS
  ss_beta <- mean(beta_ss)
  ss_gamma <- mean(gamma_ss)
  ss_alpha_sum <- alpha_ss %>%
    dplyr::select(-group) %>%
    group_by(rep) %>%
    summarise_all(sum) %>%
    dplyr::select(ss) %>%
    summarise_all(mean) %>%
    as.double()

  # SDiv
  sdiv_beta <- mean(beta_sdiv)
  sdiv_gamma <- mean(gamma_sdiv)
  sdiv_alpha_values <- alpha_sdiv %>%
    rename(sdiv_alpha = sdiv) %>%
    dplyr::select(-rep) %>%
    group_by(group) %>%
    summarise_all(mean) %>%
    ungroup() %>%
    dplyr::select(-group)
  sdiv_alpha_df <- SpatialPixelsDataFrame(xy_plots, sdiv_alpha_values, proj4string = CRS(proj4string(cube) ))
  sdiv_alpha_raster <- raster(sdiv_alpha_df)
  sdiv_alpha_mean <- mean(sdiv_alpha_values$sdiv_alpha)

  # Prepare outputs
  ss <- tibble(source = c('alpha', 'beta', 'gamma'), sum_squares = c(ss_alpha_sum, ss_beta, ss_gamma) ) %>%
    mutate(prop_gamma = sum_squares / ss_gamma)
  sdiv <- c(sdiv_alpha_mean, sdiv_beta, sdiv_gamma) ; names(sdiv) <- c('mean_alpha', 'beta', 'gamma')
  fcsd <- bind_cols(source = c('mean_alpha', 'beta', 'gamma'), bind_rows(fcsd_alpha_mean, fcsd_beta, fcsd_gamma))
  rasters <- list(beta_lcsd = lcsd_beta_raster, beta_lcss = lcss_beta_raster, alpha_sdiv = sdiv_alpha_raster, alpha_fcsd = fcsd_alpha_brick)
  out <- list(ss = ss, sdiv = sdiv, fcsd = fcsd, rasters = rasters)
  return(out)
}
