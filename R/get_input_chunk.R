#' get input chunk data from a list of rasters
#'
#' @param r_in list. list of input SpatRasters
#' @param blk list. starting row for reading and nb of rows to read
#' @param selected_bands numeric. rank of bands to select
#' @param window_size numeric. window size
#'
#' @return list. includes filtered input_data, selected_bands & nb_windows
#' @export

get_input_chunk <- function(r_in, blk, selected_bands, window_size){

  input_data <- res_shape_chunk <- list()
  nameVars <- c()
  for (fid in names(r_in)){
    input_data[[fid]] <- terra::readValues(r_in[[fid]], row = blk$row,
                                           nrows = blk$nrows, dataframe = TRUE)
    if (fid == 'mask')
      names(input_data[[fid]]) <- 'mask'
    if (dim(r_in[[fid]])[3]==1)
      nameVars <- c(nameVars, fid)
    if (dim(r_in[[fid]])[3]>1)
      nameVars <- names(input_data[[fid]])
  }
  if (is.null(selected_bands)){
    if ('mask'%in%nameVars)
      selected_bands <- seq_len(length(nameVars[-which(nameVars=='mask')]))
    if (!'mask'%in%nameVars)
      selected_bands <- seq_len(length(nameVars))
  }
  input_data <- do.call(cbind, input_data)
  names(input_data) <- nameVars
  input_data$win_ID <- produce_win_ID(inputdata = input_data, blk = blk,
                                      window_size = window_size)
  nb_windows <- max(input_data$win_ID)
  # 2a- eliminate masked pixels
  if ('mask' %in% names(input_data)){
    input_data <- input_data %>% dplyr::filter(input_data$mask > 0)
    input_data$mask <- NULL
  }
  # 2b- eliminate NA and inf
  input_data <- clean_NAsInf(input_data)
  return(list('input_data' = input_data,
              'selected_bands' = selected_bands,
              'nb_windows' = nb_windows))
}
