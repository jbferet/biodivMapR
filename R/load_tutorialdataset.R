#' loads tutorial dataset
#'
#' @param tmpdir character. path for directory
#'
#' @return Sel_PC
#' @importFrom zip unzip
#' @importFrom utils download.file
#' @export

load_tutorialdataset <- function(tmpdir = NULL){

  url_gitlab <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample'
  # url for a S2 subset
  name_raster <- 'S2A_T33NUD_20180104_Subset'
  urlraster <- file.path(url_gitlab,'RASTER',name_raster)
  # create a temporary directory (choose your own data directory)
  if (is.null(tmpdir)) tmpdir <- 'biodivMapR_tutorial_data'
  dir.create(path = tmpdir, showWarnings = F, recursive = T)
  # create destination raster file
  path_raster <- file.path(tmpdir,name_raster)
  download.file(url = urlraster, destfile = path_raster,
                method = 'auto', quiet = FALSE, mode = "wb")

  # url for the S2 subset header
  urlhdr <- paste0(urlraster,'.hdr')
  # name your raster HDR with the same name as the binary raster, with .hdr extension
  path_hdr <- paste0(path_raster,'.hdr')
  download.file(url = urlhdr, destfile = path_hdr,
                method = 'auto', quiet = FALSE, mode = "w")

  # url for vector data
  name_zip_vect <- 'S2A_T33NUD_Plots.zip'
  urlvect <- file.path(url_gitlab,'VECTOR',name_zip_vect)
  # name zip file including plots located on the tile
  path_zip <- file.path(tmpdir,name_zip_vect)
  download.file(url = urlvect, destfile = path_zip)
  destunz <- file.path(tmpdir,tools::file_path_sans_ext(name_zip_vect))
  zip::unzip(zipfile = path_zip,exdir = destunz)
  return(list('path_raster' = path_raster,
              'path_hdr' = path_hdr,
              'path_vector' = destunz))
}
