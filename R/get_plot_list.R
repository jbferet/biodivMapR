#' get a list of plots defined in a vector file defined as input
#'
#' @param dsn character. path for vector file
#' @param nbdigits numeric. number of digits to define plot names
#'
#' @return list of polygons
#' @importFrom sf st_read
#' @export

get_plot_list <- function(dsn, nbdigits = 3){
  plots_sf <- sf::st_read(dsn = dsn, quiet = T)[[1]]
  nbPlots <- length(plots_sf)
  plots <- list()
  for (i in seq_len(nbPlots))
    plots[[i]] <- plots_sf[i]
  names(plots) <- num2char(seq_len(nbPlots),
                           nbdigits = nbdigits)
  return(plots)
}
