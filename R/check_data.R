#' Checks if the data to be processed has the format type expected
#'
#' @param input_data character. which check?
#' @param arguments list. additional arguments
#'
#' @return None
#' @importFrom terra names
#' @export

check_data <- function(input_data, arguments) {

  if (inherits(arguments,'character')){
    if (arguments=='input_rast'){
      wl <- as.numeric(terra::names(input_data))
      if (length(which(is.na(wl)))>0){
        whichNA <- which(is.na(as.numeric(terra::names(input_data))))
        message("***************************************************************")
        message("     No wavelength is associated to SpatRaster 'input_rast'.   ")
        message(" Please specify central wavelength (in nm) as names(input_rast)")
        message("         if processing multi / hyperspectral optical data      ")
        message("   Otherwise, set Continuum_Removal = FALSE and do not use     ")
        message("                   radiometric filtering                       ")
        message("***************************************************************")
        print(terra::names(input_data)[whichNA])
      }
    }
    if (arguments=='input_mask'){
    }
  } else if (inherits(arguments,'list')){
    if (!is.null(arguments$filter)){
      if (arguments$filter == 'WL'){
        # set wl
        wl <- as.numeric(input_data)
        # Dist2Band = distance to theoretical band for blue, red and NIR data
        if (max(wl)<40)
          stop('Define wavelength for input_rast in nanometers')
      }
    }
  }
  return(invisible())
}
