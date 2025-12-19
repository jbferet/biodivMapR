#' Check if principal components are properly selected as expected by the method
#'
#' @param pca_rast_path character. Path of the PCA image
#' @param output_dir character. Path for output directory
#' @param File_Open Boolean. Set to TRUE for file to open automatically
#'
#' @return Sel_PC
#' @importFrom utils file.edit read.csv
#' @export

select_PCA_components <- function(pca_rast_path,
                                  output_dir,
                                  File_Open = TRUE) {

  message("Please check following PCA file:")
  print(pca_rast_path)
  message("and identify the components selected to compute spectral diversity")
  message("one component per line: carriage return after each component")
  Sel_PC <- file.path(output_dir, 'Selected_Components.txt')
  message("")
  message("Then press 'Save'")
  if (!file.exists(Sel_PC))
    file.create(Sel_PC)
  if (File_Open == TRUE)
    utils::file.edit(Sel_PC, title=basename(Sel_PC))
  message("")
  message("Selected components were saved here")
  print(Sel_PC)
  message("Please press Enter to continue the process")
  readline(prompt = "")
  Selected_PC <- utils::read.csv(file = Sel_PC, header = FALSE)[[1]]
  return(Selected_PC)
}
