#' remove constant bands
#'
#' @param DataMatrix numeric. each variable is a column
#' @param Spectral list. summary of spectral information: which spectral bands selected from initial data
#'
#' @return updated DataMatrix and Spectral
#' @importFrom stats sd
#' @export

rm_invariant_bands <- function(DataMatrix, Spectral) {
  # samples with inf value are eliminated
  for (i in 1:ncol(DataMatrix)) {
    elim <- which(DataMatrix[, i] == Inf)
    if (length(elim) > 0) DataMatrix <- DataMatrix[-elim, ]
  }
  # bands showing null std are removed
  stdsub <- apply(DataMatrix, 2, sd)
  BandsNoVar <- which(stdsub == 0)
  # BandsNoVar  = which(stdsub<=0.002)
  if (!length(Spectral$Bands2Keep[BandsNoVar]) == 0) {
    DataMatrix <- DataMatrix[, -BandsNoVar]
  }
  # !! the wl which is discarded correspond to original spectral bands,
  # whereas BandsNoVar corresponds to spectral band after contiuum removal
  Spectral$BandsNoVar <- BandsNoVar
  return(list("DataMatrix" = DataMatrix, "Spectral" = Spectral))
}
