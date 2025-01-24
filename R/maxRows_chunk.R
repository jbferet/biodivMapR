#' redefined chunks based on the max number of rows per chunk
#
#' @param blk list. obtained from function terra::blocks
#' @param maxRows numeric. max number of rows in each block
#
#' @return blk
#' @export

maxRows_chunk <- function(blk, maxRows = NULL){
  # adjust number of lines to be read each time
  if (!is.null(maxRows)){
    if (blk$nrows[1]>maxRows){
      nrows_indiv <- maxRows
      blkud <- list()
      blkud$nrows <- rep(nrows_indiv,floor(sum(blk$nrows)/nrows_indiv))
      if (sum(blkud$nrows)<sum(blk$nrows)){
        blkud$nrows <- c(blkud$nrows, sum(blk$nrows)-sum(blkud$nrows))
      }
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

