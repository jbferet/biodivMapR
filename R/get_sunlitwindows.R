#' eliminate windows with insufficient sunlit pixels
#
#' @param inputdata list. input data with window_ID produced from produce_win_ID
#' @param pixperplot numeric. minimum number of pixels per window
#' @param MinSun numeric. minimum percentage of sunlit pixels
#
#' @return NMDS position based on nearest neighbors from NMDS
#' @importFrom dplyr group_by relocate %>% last_col
#' @importFrom purrr map
#' @importFrom tidyr unnest
#' @export
#'
get_sunlitwindows <- function(inputdata, pixperplot, MinSun = 0.25){
  inputwindow <- inputdata %>% dplyr::group_by(win_ID) %>% nest()
  nbPix_Sunlit <- unlist(purrr::map(inputwindow$data,nrow))
  PCsun <- nbPix_Sunlit/pixperplot
  SelWindows <- which(PCsun > MinSun)
  inputwindow <- inputwindow[SelWindows,]
  inputwindow <- inputwindow %>% unnest(win_ID) %>% unnest(data)
  inputwindow <- inputwindow %>% relocate(win_ID, .after = last_col())
  return(inputwindow)
}
