terra_extract <- function(input_rast, y, p = NULL){
  rast_sample <- terra::extract(x = input_rast, y = y)
  if (!is.null(p)) p()
  return(rast_sample)
}
