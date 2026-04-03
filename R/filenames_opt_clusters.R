#' checks if files produced from biodivMapR_opt_cluster exist
#'
#' @param outputdir character
#' @param obs_criterion character
#'
#' @return files_exist boolean if all expected files exist
#'
#' @export
filenames_opt_clusters <- function(outputdir, obs_criterion){

  file_exists <- c()
  for (criterion in obs_criterion){
    filename_pearson_mean <- file.path(outputdir, paste0(criterion, '_pearson_mean.csv'))
    filename_spearman_mean <- file.path(outputdir, paste0(criterion, '_spearman_mean.csv'))
    filename_pearson_all <- file.path(outputdir, paste0(criterion, '_pearson_all.csv'))
    filename_spearman_all <- file.path(outputdir, paste0(criterion, '_spearman_all.csv'))
    if (file.exists(filename_pearson_mean) & file.exists(filename_spearman_mean) &
        file.exists(filename_pearson_all) & file.exists(filename_spearman_all)){
      file_exists <- c(file_exists, TRUE)
    } else {
      file_exists <- c(file_exists, FALSE)
    }
  }
  files_exist <- all(file_exists)
  return(files_exist)
}
