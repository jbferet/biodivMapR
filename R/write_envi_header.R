#' writes ENVI hdr file
#'
#' @param hdr content to be written
#' @param hdr_path Path of the hdr file
#'
#' @return None
#' @export
write_envi_header <- function(hdr, hdr_path) {
  h <- lapply(hdr, function(x) {
    if (length(x)>1 || (is.character(x) && stringr::str_count(x, "\\w+") > 1)){
      x <- paste0("{", paste(x, collapse = ","), "}")
    }
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c("ENVI", paste(names(hdr), h, sep = " = ")), con = hdr_path)
  return(invisible())
}
