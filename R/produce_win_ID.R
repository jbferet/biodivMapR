#' design a matrix with window ID based on an original raster and window size in pixels
#'
#' @param inputdata numeric. image chunk
#' @param blk list.
#' @param window_size numeric. window size for square plots
#'
#' @return blk
#' @export

produce_win_ID <- function(inputdata, blk, window_size){
  # generate win_ID crresponding to an image, given chunk & window size
  nbRows <- ceiling(blk$nrows/window_size)
  nbRows_piece <- rep(x = window_size,nbRows)
  if (blk$nrows%%window_size > 0) nbRows_piece[nbRows] <- blk$nrows%%window_size
  nbCols <- ceiling((nrow(inputdata)/blk$nrows)/window_size)
  win_ID <- win_ID_lowRes <- list()
  for (i in 1:nbRows){
    # get window ID for 1 line of future product at original res x window_size res
    win_ID_lowRes[[i]] <- seq(1,nbCols)
    # get window ID for 1 line of current raster at original res
    widthwin_ID <- matrix(rep(win_ID_lowRes[[i]],
                              each = window_size)[1:(nrow(inputdata)/blk$nrows)],
                          nrow = 1)
    # get window ID for a chunk of current raster at original res
    Chunkwin_ID <- widthwin_ID %x% rep(1, nbRows_piece[i])
    win_ID[[i]] <-  c(t(Chunkwin_ID))+(i-1)*nbCols
  }
  win_ID <- unlist(win_ID)
  win_ID_lowRes <- unlist(win_ID_lowRes)
  return(win_ID)
}

