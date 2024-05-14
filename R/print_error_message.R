#' prints an error message if problems occur
#'
#' @param def_error character. nature of the error
#' @param optarg character. optional additional info required for the error message
#'
#' @return none.
#' @export

print_error_message <- function(def_error, optarg) {
  if (def_error=='error_input'){
    message("")
    message("*********************************************************")
    message("WARNING: the processing resulted in NA or infinite values")
    message("     This may be due to noisy spectral domains or        ")
    message(" individual pixels showing Inf or Na values in input data")
    message("               Please check input data                   ")
    message("                                                         ")
    message("   if nothing wrong identified with input data, please   ")
    message("   submit a bug report, reproduce the bug with reduced   ")
    message("     dataset and contact the authors of the package      ")
    message("                   process aborted                       ")
    message("*********************************************************")
    message("")
  }
  if (def_error=='error_no_PCA_file'){
    PCA_Files <- optarg
    message("")
    message("*********************************************************")
    message("WARNING: This file required to compute spectral species is missing")
    print(PCA_Files)
    message("process aborted")
    message("*********************************************************")
    message("")
    stop()
  }
  if (def_error=='error_PC_sel'){
    Output_Dir_PCA <- optarg
    print("PC SELECTION MUST BE PERFORMED FIRST")
    print("Please identify selected components either in this file:")
    PC_Select_Path <- file.path(Output_Dir_PCA, "Selected_Components.txt")
    print(PC_Select_Path)
    print("or in the 'SelectedPCs' variable of map_spectral_species")
    print("Image processing aborted")
    stop()
  }
  if (def_error=='missing_hdr'){
    ImPathHDR <- optarg
    message("WARNING : COULD NOT FIND HDR FILE ", ImPathHDR)
  }
  if (def_error=='Beta_info_file_missing'){
    message('Warning init_PCoA: no file corresponding to input variable "Beta_info_read"')
    message('the initialization for beta diversity mapping will be performed')
  }
  if (def_error=='Kmeans_info_file_missing'){
    message('Warning init_Kmeans: no file corresponding to input variable "Kmeans_info_read"')
    message('the initialization for Kmeans clustering will be performed')
  }
  return(invisible())
}
