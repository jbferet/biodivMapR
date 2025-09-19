#' Computes continuum removal for matrix shaped data: more efficient than
#' processing individual spectra
#' the convex hull is based on the computation of the derivative between R at a
#' given spectral band and R at the following bands
#'
#' @param mat_init numeric. initial data matrix (nb samples x nb bands)
#' @param spectral_bands numeric. central wavelength for the spectral bands
#' @param p list. progressor object for progress bar
#
#' @return samples from image and updated number of pixels to sampel if necessary
#' @export

continuum_removal <- function(mat_init, spectral_bands, p = NULL) {

  # Filter and prepare data prior to continuum removal
  cr_data <- filter_prior_cr(mat_init, spectral_bands)
  nb_bands <- dim(cr_data$mat_init)[2]
  nb_samples <- cr_data$nb_samples
  # mat_init <- cr_data$mat_init
  # cr_data$mat_init <- c()
  # spectral_bands <- cr_data$spectral_bands
  nb_samples_update <- length(cr_data$samples_to_keep)

  # if samples to be considered
  if (nb_samples > 0) {
    cr_results <- matrix(0, ncol = (nb_bands - 2), nrow = nb_samples)
    cr_results0 <- continuumRemoval(X = cr_data$mat_init,
                                    wav = cr_data$spectral_bands)
    cr_results[cr_data$samples_to_keep, ] <- as.matrix(cr_results0[,2:(nb_bands-1)])
    # # initialization:
    # # spectral band corresponding to each element of the data matrix
    # lambda <- repmat(matrix(spectral_bands, nrow = 1), nb_samples_update, 1)
    # # prepare matrices used to check evolution of the CR process
    # # - elements still not processed through continuum removal: initialization to 1
    # still_need_cr <- matrix(1, nrow = nb_samples_update, ncol = nb_bands)
    # # - value of the convex hull: initially set to 0
    # convex_hull <- matrix(0, nrow = nb_samples_update, ncol = nb_bands)
    # # - reflectance value for latest interception with convex hull:
    # # initialization to value of the first reflectance measurement
    # intercept_hull <- repmat(matrix(mat_init[, 1], ncol = 1), 1, nb_bands)
    # # - spectral band of latest interception
    # latest_intercept <- repmat(X = matrix(spectral_bands[1], ncol = 1),
    #                            m = nb_samples_update, n = nb_bands)
    # # number of spectral bands found as intercept
    # nb_intercept <- 0
    # # continues until arbitrary stopping criterion:
    # # stops when reach last spectral band (all values before last = 0)
    # # while (max(still_need_cr[, seq_len(nb_bands - 2)]) == 1 & (nb_intercept <= (nb_bands / 2))) {
    # while (max(still_need_cr[, seq_len((nb_bands - 2))]) == 1) {
    #   nb_intercept <- nb_intercept + 1
    #   # identify samples still needing continuum removal
    #   sel <- which(still_need_cr[,(nb_bands-2)]==1)
    #   # update variables to process samples needing CR only
    #   nb_samples_update_tmp <- length(sel)
    #   lambda_tmp <- lambda[sel,]
    #   mat_init_tmp <- mat_init[sel,]
    #   latest_intercept_tmp <- latest_intercept[sel,]
    #   still_need_cr_tmp <- still_need_cr[sel,]
    #   convex_hull_tmp <- convex_hull[sel,]
    #   intercept_hull_tmp <- intercept_hull[sel,]
    #   # Mstep give the position of the values to be updated
    #   update_data <- matrix(1, nrow = nb_samples_update_tmp, ncol = nb_bands)
    #   update_data[, nb_bands] <- 0
    #   # initial step: first column set to 0; following steps: all bands below
    #   # max of the convex hull are set to 0
    #   update_data[which((lambda_tmp - latest_intercept_tmp) < 0)] <- 0
    #   # compute slope for each coordinate
    #   slope <- as.matrix((mat_init_tmp - intercept_hull_tmp) / (lambda_tmp - latest_intercept_tmp) * still_need_cr_tmp)
    #   # set current spectral band and previous bands to -9999
    #   if (!length(which(still_need_cr_tmp == 0)) == 0) {
    #     slope[which(still_need_cr_tmp == 0)] <- -9999
    #   }
    #   if (!length(which(is.na(slope))) == 0) {
    #     slope[which(is.na(slope))] <- -9999
    #   }
    #   # get max index for each row and convert into linear index
    #   index_max_slope <- RowToLinear(max.col(slope, ties.method = "last"),
    #                                  nb_samples_update_tmp, nb_bands)
    #   # !!!! OPTIM: replace repmat with column operation
    #   # update coordinates of latest intercept
    #   latest_intercept_tmp <- repmat(matrix(lambda_tmp[index_max_slope],
    #                                         ncol = 1), 1, nb_bands)
    #   # update latest intercept
    #   intercept_hull_tmp <- repmat(matrix(as.matrix(mat_init_tmp)[index_max_slope],
    #                                       ncol = 1), 1, nb_bands)
    #   # values corresponding to the domain between the two continuum maxima
    #   update_data[which((lambda_tmp - latest_intercept_tmp) >= 0 |
    #                       latest_intercept_tmp == spectral_bands[nb_bands])] <- 0
    #   # values to eliminate for the next analysis: all spectral bands before latest intercept
    #   still_need_cr_tmp[which((lambda_tmp - latest_intercept_tmp) < 0)] <- 0
    #   # the max slope is known, as well as the coordinates of the beginning and ending
    #   # a matrix now has to be built
    #   convex_hull_tmp <- convex_hull_tmp +
    #     update_data * (intercept_hull_tmp + sweep((lambda_tmp - latest_intercept_tmp),
    #                                               1, slope[index_max_slope], "*"))
    #   # update variables
    #   convex_hull[sel,] <- convex_hull_tmp
    #   still_need_cr[sel,] <- still_need_cr_tmp
    #   lambda[sel,] <- lambda_tmp
    #   latest_intercept[sel,] <- latest_intercept_tmp
    #   intercept_hull[sel,] <- intercept_hull_tmp
    # }
    # cr_results0 <- mat_init[, 2:(nb_bands - 2)] / convex_hull[, 2:(nb_bands-2)]
    # cr_results <- matrix(0, ncol = (nb_bands - 3), nrow = nb_samples)
    # cr_results[cr_data$samples_to_keep, ] <- as.matrix(cr_results0)
  } else {
    cr_results <- matrix(0, ncol = (nb_bands - 3), nrow = nb_samples)
  }
  cr_results <- data.frame(cr_results)
  if (!is.null(p))
    p()
  list <- ls()
  rm(list = list[-which(list == "cr_results")])
  rm(list)
  gc()
  return(cr_results)
}
