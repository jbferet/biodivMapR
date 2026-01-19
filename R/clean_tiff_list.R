#' computes diversity metrics from raster
#'
#' @param list_files vector. list of file names corresponding to tiff
#'
#' @return list_files vector. list of file names withoout aux.xml and tfw files
#' @export

clean_tiff_list <- function(list_files){
  list_files <- unique(gsub(pattern = '.aux.xml', replacement = '',
                            x = list_files))
  list_files <- unique(gsub(pattern = '.tfw', replacement = '.tiff',
                            x = list_files))
  return(list_files)
}
