#' This function gets path from an asset in the JSON file
#'
#' @param jsonobj character. path for json file
#' @param asset character. name of the asset of interest
#' @param dateAcq date. date of acquisition
#'
#' @return asset_path
#' @export

get_asset_path <- function(jsonobj, asset, dateAcq){
  listAssets <- unlist(lapply(lapply(lapply(jsonobj$features,'[[','assets'),
                                     '[[',asset),'[[','href'))
  DatesAcq <- get_date(prodName = basename(listAssets))
  seldate <- which(DatesAcq==dateAcq)
  asset_path <- listAssets[seldate]
  return(asset_path)
}
