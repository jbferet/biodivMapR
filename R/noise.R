#' this function
#'
#' @param X numeric. data matrix
#' @param coordPix numeric.
#'
#' @return rescaled data, min and max values
#' @export
# If coordPix is not NULL, X and coordPix are exepected to have the same order,
# i.e. coordPix[1, ] corresponds to X[1, ], coordPix[2, ] corresponds to X[2, ], ...
noise <- function(X, coordPix=NULL){
  if(is.null(coordPix)){
    # if(matlab::ndims(X)!=3)
    if(length(dim(X))!=3)
      stop('X is expected to be a 3D array: y,x,band for row,col,depth.')
    Xdim <- dim(X)
    # Shift x/y difference
    Y = ((X[2:Xdim[1],2:Xdim[2],]-X[1:(Xdim[1]-1),2:Xdim[2],]) +
           (X[2:Xdim[1],2:Xdim[2],]-X[2:Xdim[1],1:(Xdim[2]-1),]))/2
  }else{
    if(!all(c('Kind', 'id') %in% colnames(coordPix)))
      stop("Columns 'Kind' and 'id' are missing in coordPix.")
    kernel = matrix(0, 3, 3)
    kernel[c(5, 6, 8)]=c(1, -1/2, -1/2)

    if(!identical(order(coordPix$id), 1:nrow(coordPix)))
      stop("coordPix is not ordered along column 'id'. Order coordPix as well as X before trying again.")

    Y=0
    for(ik in which(kernel!=0)){
      Y = Y + X[coordPix$Kind==ik,]*kernel[ik]
    }
  }
  return(Y)
}
