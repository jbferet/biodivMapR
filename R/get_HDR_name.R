#' get hdr name from image file name
#'
#' @param ImPath character. ath of the image
#' @param showWarnings boolean. set TRUE if warning because HDR does not exist
#'
#' @return corresponding hdr
#' @export

get_HDR_name <- function(ImPath,showWarnings=TRUE) {
  ImPathHDR <- paste0(file_path_sans_ext(ImPath), ".hdr")
  if (showWarnings==TRUE){
    if (!file.exists(ImPathHDR))
      print_error_message(def_error = 'missing_hdr', optarg = ImPathHDR)
  }
  return(ImPathHDR)
}
