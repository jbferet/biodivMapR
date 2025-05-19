#' adjusts number of rows from chunks
#'
#' @param blk list.
#' @param nbRows numeric
#'
#' @return blk
#' @export

nbRows_chunk <- function(blk, nbRows = NULL){
  # adjust number of lines to be read each time
  if (!is.null(nbRows)){
    if (length(blk$nrows) >1){
      nrows_indiv <- blk$nrows[1] - (blk$nrows[1] %% nbRows)
      if (nrows_indiv==0)
        nrows_indiv <- nbRows
      blkud <- list()
      blkud$nrows <- rep(nrows_indiv,floor(sum(blk$nrows)/nrows_indiv))
      if (sum(blkud$nrows)<sum(blk$nrows))
        blkud$nrows <- c(blkud$nrows, sum(blk$nrows)-sum(blkud$nrows))
      blkud$row <- NULL
      blkud$n <- length(blkud$nrows)
      blkud$row <- c(1,cumsum(blkud$nrows)+1)[seq_len(blkud$n)]
      blk$row <- blkud$row
      blk$nrows <- blkud$nrows
      blk$n <- blkud$n
    }
  }
  return(blk)
}

