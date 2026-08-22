#' Harmonize ENVI header file
#'
#' @param output_raster character. path for raster file
#' @param filetype character. raster driver
#' @param bandNames character. name for bands in raster
#'
#' @return none
#' @export

harmonize_envi_hdr <- function(output_raster, filetype, bandNames){
  if (filetype == 'EHdr')
    raster::hdr(raster::raster(output_raster), format = "ENVI")
  hdr <- read_envi_header(get_hdr_name(image_path = output_raster))
  if (filetype == 'EHdr')
    hdr$`coordinate system string` <- utils::read.table(paste0(file_path_sans_ext(output_raster), ".prj"))
  if (!is.null(bandNames) &
      dim(raster::stack(output_raster))[3] == length(bandNames)) {
    hdr$`band names` <- paste0('',bandNames,'')
  }
  write_envi_header(hdr = hdr, hdr_path = get_hdr_name(output_raster))
}
