#' converts numeric value to characters with specific lengths, adding 0 at beginning
#'
#' @param val numeric. value to be converted into string
#' @param nbdigits numeric. number pf digits for string
#'
#' @return valchar string corresponding to val with 0 if needed
#' @export
#'
num2char <- function(val, nbdigits = 3){
  valchar <- as.character(val)
  for (i in seq_len(length(valchar))){
    if (nchar(valchar[i])<nbdigits){
      fact <- 10**seq(from = nchar(valchar[i]), to = (nbdigits-1))
      for (fact10 in fact) valchar[i] <- paste0('0',valchar[i])
    }
  }
  return(valchar)
}
