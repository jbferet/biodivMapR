#' sample pixels from raster data
#'
#' @param input_rast SpatRaster or list of SpatRaster
#' @param pix2extract dataframe. lat and lon for pixels to extract
#' @param xy boolean. set TRUE to get coordinates for each pixel sampled
#' @param prog boolean. set TRUE to display progression bar
#'
#'
#' @return rast_sample dataframe. pixel/plot info extracted from input_rast
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom terra extract
#' @importFrom stats na.omit
#' @export

sample_raster <- function(input_rast, pix2extract, xy = FALSE, prog = FALSE){
  # if multilayer spatial raster provided as input
  if (inherits(x = input_rast,what = 'SpatRaster')){
    if (inherits(pix2extract,what = 'matrix')){
      rast_sample <- terra::extract(x = input_rast, y = pix2extract)
    } else {
      rast_sample <- terra::extract(x = input_rast, y = pix2extract, xy = xy)
    }
    # if list of rasters
  } else if (inherits(x = input_rast, what = 'list')){
    rast_sample <- NULL
    names_layers <- names(input_rast)
    if (prog){
      # progressr::handlers(global = TRUE)
      suppressWarnings(progressr::handlers("cli"))
      # progressr::handlers("debug")
      suppressWarnings(with_progress({
        p <- progressr::progressor(steps = length(input_rast))
        rast_sample_list <- lapply(X = input_rast,
                                   FUN = terra_extract,
                                   y = pix2extract,
                                   p = p)
      }))
    } else {
      rast_sample_list <- lapply(X = input_rast,
                                 FUN = terra_extract,
                                 y = pix2extract)
    }
    dimsextract <- lapply(rast_sample_list,ncol)
    selcol <- 1
    if (inherits(pix2extract,what = 'SpatVector')) selcol <- 2
    if (all(unlist(dimsextract)==selcol)){
      rast_sample <- data.frame(lapply(rast_sample_list, '[[',selcol))
    } else {
      sel1 <- which(unlist(dimsextract)==selcol)
      rast_sample <- data.frame(lapply(rast_sample_list[sel1], '[[',selcol))
      sel <- which(unlist(dimsextract)>selcol)
      for (mat in sel){
        if (dim(rast_sample)[1]>0){
          if (inherits(rast_sample_list[sel],what = 'list')){
            rast_sample <- cbind(rast_sample,rast_sample_list[sel][[1]])
            rast_sample$ID <- NULL
          } else if (inherits(rast_sample_list[sel],what = 'data.frame')){
            rast_sample <- cbind(rast_sample,rast_sample_list[sel])
            rast_sample$ID <- NULL
          }
        } else {
          if (inherits(rast_sample_list[sel],what = 'list')){
            rast_sample <- rast_sample_list[sel][[1]]
          } else if (inherits(rast_sample_list[sel],what = 'data.frame')){
            rast_sample <- rast_sample_list[sel]
          }
        }
      }
    }

    if (!is.null(rast_sample_list[[1]]$ID)) rast_sample$ID <- rast_sample_list[[1]]$ID
    #   # if list of individual layers
    #   names_before <- names(rast_sample)
    #   if (dim(input_rast[[lr]])[3]==1) {
    #     rast_sample <- cbind(rast_sample, rast_sample0[,2])
    #     names(rast_sample) <- c(names_before, names_layers[lr])
    #   } else if (dim(input_rast[[lr]])[3]>1) {
    #     rast_sample <- cbind(rast_sample, rast_sample0[,-1])
    #     names(rast_sample) <- c(names_before, names(input_rast[[lr]]))
    #   }
    # for (lr in seq_len(length(input_rast))){
    #   rast_sample0 <- terra::extract(x = input_rast[[lr]], y = pix2extract)
    #   if (is.null(rast_sample)) rast_sample <- data.frame('ID' = rast_sample0$ID)
    #   # if list of individual layers
    #   names_before <- names(rast_sample)
    #   if (dim(input_rast[[lr]])[3]==1) {
    #     rast_sample <- cbind(rast_sample, rast_sample0[,2])
    #     names(rast_sample) <- c(names_before, names_layers[lr])
    #   } else if (dim(input_rast[[lr]])[3]>1) {
    #     rast_sample <- cbind(rast_sample, rast_sample0[,-1])
    #     names(rast_sample) <- c(names_before, names(input_rast[[lr]]))
    #   }
    #   pb$tick()
    # }
  }
  return(rast_sample)
}

