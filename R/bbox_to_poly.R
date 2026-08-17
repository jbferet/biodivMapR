#' converts a bbox to a polygon
#'
#' @param x bbox
#' @param crs crs of the polygon
#'
#' @return polyg polygon
#' @importFrom sf st_polygon st_sfc
#' @export

bbox_to_poly <- function(x, crs = 4326) {
  xmin <- x[["xmin"]]
  ymin <- x[["ymin"]]
  xmax <- x[["xmax"]]
  ymax <- x[["ymax"]]
  polyg <- rbind(c(xmin, ymin), c(xmax, ymin),
        c(xmax, ymax), c(xmin, ymax),
        c(xmin, ymin)) |>
    list() |>
    sf::st_polygon() |>
    sf::st_sfc(crs = crs)
  return(polyg)
}
