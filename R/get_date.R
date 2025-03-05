#' gets acquisition date from S2 image
#'
#' @param prodName character. original name for the S2 image
#'
#' @return DateAcq character
#' @export

get_date <- function(prodName){
  prodName <- basename(prodName)
  sel1 <- gsub(pattern = 'R\\d{3}', replacement = '', x = prodName)
  sel2 <- gsub(pattern = '.*MSIL2A_', replacement = '', x = sel1)
  DateAcq <- as.Date(gsub('T.*', '', sel2),format = '%Y%m%d')
  return(DateAcq)
}
