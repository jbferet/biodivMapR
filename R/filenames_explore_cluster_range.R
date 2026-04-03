#' checks if files produced from explore_cluster_range exist
#'
#' @param outputdir character
#' @param obs_criterion character
#' @param repet character
#'
#' @return files_exist boolean if all expected files exist
#'
#' @export
filenames_explore_cluster_range <- function(outputdir, obs_criterion, repet){

  file_exists <- c()
  for (criterion in obs_criterion){
    filename <- file.path(outputdir, paste0(criterion,'_',repet,'.csv'))
    if (file.exists(filename)){
      file_exists <- c(file_exists, TRUE)
    } else {
      file_exists <- c(file_exists, FALSE)
    }
  }
  files_exist <- all(file_exists)
  return(files_exist)
}
