#' computes k-means from nbIter subsets taken from dataPCA
#'
#' @param rast_sample numeric. initial dataset sampled from PCA image
#' @param nbIter numeric. nb of iterations averaged to compute diversity indices
#' @param nbclusters numeric. number of clusters used in kmeans
#' @param nbCPU numeric. Number of CPUs available
#' @param algorithm character. algorithm used in the kmeans clustering
#' @param progressbar boolean. set true for progress bar during clustering
#'
#' @return list of centroids and parameters needed to center/reduce data
#' @import cli
#' @importFrom progressr progressor handlers with_progress
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats kmeans
#' @importFrom snow splitRows
#' @importFrom parallel makeCluster stopCluster
#'
#' @export

get_kmeans <- function(rast_sample, nbIter, nbclusters = 50,
                       nbCPU = 1, algorithm = 'Hartigan-Wong',
                       progressbar = TRUE) {
  # define boundaries defining outliers based on IQR
  IQRsubset <- lapply(X = rast_sample, IQR_outliers, weightIRQ = 2)
  m0 <- unlist(lapply(IQRsubset,'[[',1))
  M0 <- unlist(lapply(IQRsubset,'[[',2))
  d0 <- M0 - m0
  if (any(is.na(c(m0, M0, d0))) | any(is.infinite(c(m0, M0, d0)))){
    print_error_message('error_input')
    return(list("Centroids" = NULL,
                "MinVal" = m0, "MaxVal" = M0, "Range" = d0,
                "Error" = TRUE))
  } else {
    rast_sample <- center_reduce(rast_sample, m0, d0)
    rast_sample <- snow::splitRows(x = rast_sample, ncl = nbIter)
    if (nbCPU>1){
      # plan(multisession, workers = nbCPU) ## Parallelize using four cores
	  cl <- parallel::makeCluster(nbCPU)
      plan("cluster", workers = cl)  ## same as plan(multisession, workers = nbCPU)

      fun_apply <- future_lapply
    } else {
      fun_apply <- lapply
    }
    if (progressbar==TRUE){
      handlers(global = TRUE)
      handlers("cli")
      with_progress({
        p <- progressr::progressor(steps = nbIter)
        res <- fun_apply(X = rast_sample,
                         FUN = kmeans_progressr,
                         centers = nbclusters,
                         iter.max = 50, nstart = 10,
                         algorithm = algorithm, p = p)
      })
    } else {
      res <- fun_apply(X = rast_sample,
                       FUN = kmeans_progressr,
                       centers = nbclusters,
                       iter.max = 50, nstart = 10,
                       algorithm = algorithm, p = NULL)
    }
    if (nbCPU>1) parallel::stopCluster(cl)
    if (nbCPU>1) plan(sequential)
    Centroids <- lapply(res,'[[',2)
    return(list("Centroids" = Centroids,
                "MinVal" = m0, "MaxVal" = M0, "Range" = d0, "Error" = FALSE))
  }
}
