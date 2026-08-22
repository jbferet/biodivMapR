#' get hdr name from image file name
#'
#' @param image_path character. ath of the image
#' @param showWarnings boolean. set TRUE if warning because HDR does not exist
#'
#' @return path for corresponding hdr
#' @export

get_hdr_name <- function(image_path,showWarnings=TRUE) {
  image_path_hdr <- paste0(file_path_sans_ext(image_path), ".hdr")
  if (showWarnings==TRUE){
    if (!file.exists(image_path_hdr))
      print_error_message(def_error = 'missing_hdr', optarg = image_path_hdr)
  }
  return(image_path_hdr)
}
