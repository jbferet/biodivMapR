#' split chunk into subsets to prepare for parallel processing
#
#' @param SSchunk  list. obtained from get_spectralSpecies
#' @param nbCPU numeric. Number of CPUs available
#
#' @return estimated NMDS position based on nearest neighbors from NMDS
#' @importFrom dplyr group_by
#' @importFrom tidyr nest
#' @export
#'
split_chunk <- function(SSchunk, nbCPU){

  win_ID <- NULL
  nbWindows <- length(unique(SSchunk$win_ID))
  windowperCPU <- ceiling(nbWindows/nbCPU)
  lbWin <- seq(1,nbWindows, by= windowperCPU)
  ubWin <- c(lbWin[-1]-1,nbWindows)
  SSwindow <- SSchunk %>% dplyr::group_by(win_ID) %>% nest()

  # warning, output is not sorted as expected!!
  sortUpdate <- sort(SSwindow$win_ID,index.return=T)
  SSwindow$data <- SSwindow$data[sortUpdate$ix]
  SSwindow$win_ID <- SSwindow$win_ID[sortUpdate$ix]

  # initialize
  SSwindow_perCPU <- IDwindow_perCPU <- list()
  for (i in seq_len(length(lbWin))) {
    SSwindow_perCPU[[i]] <- SSwindow$data[lbWin[i]:ubWin[i]]
    IDwindow_perCPU[[i]] <- SSwindow$win_ID[lbWin[i]:ubWin[i]]
  }
  return(list('SSwindow_perCPU' = SSwindow_perCPU,
              'IDwindow_perCPU' = IDwindow_perCPU))
}
