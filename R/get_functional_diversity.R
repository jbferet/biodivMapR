#' get functional diversity metrics from dataframe
#' This function was inspired from FD package
#' @param spectraits numeric. dataframe containing species in rows and trait values in columns
#' @param FDmetric character. Functional diversity metric
#' @param p list. progressor object for progress bar
#
#' @return FDmetrics
#' @importFrom geometry convhulln
#' @importFrom ape mst
#' @importFrom stats dist
#'
#' @export

get_functional_diversity <- function(spectraits,
                                     FDmetric = c('FRic', 'FEve', 'FDiv'),
                                     p = NULL){
  nbTraits <- ncol(spectraits)
  nbSpecies <- nrow(spectraits)
  FRic <- FDiv <- FEve <- NA

  # if (FDmetric %in% c('FRic', 'FDiv')){
  if ('FRic' %in% FDmetric | 'FDiv' %in% FDmetric){
    if (nbTraits==1){
      FRic <- max(spectraits, na.rm = T) - min(spectraits, na.rm = T)
      FDiv <- sum((spectraits[[1]]-mean(spectraits[[1]],na.rm = T))**2, na.rm = T)**0.5
    }
    if (nbTraits>1){
      # if (!is.na(match('FRic', FDmetric)) | !is.na(match('FDiv', FDmetric))){
      # inspired from FD package
      # convex hull using geometry
      FunctHull <- geometry::convhulln(spectraits, "Fx TO 'vert.txt'", output.options = 'FA')
      # 1- Functional Richness
      if (!is.na(match('FRic', FDmetric))){
        FRic <- FunctHull$vol
      }
      # 2- Functional Divergence mean distance from centroid
      if (!is.na(match('FDiv', FDmetric))){
        vert1 <- scan("vert.txt", quiet = T)
        vert2 <- vert1 + 1
        vertices <- vert2[-1]
        trvertices <- spectraits[vertices, ]
        # coordinates of the center of gravity of the vertices (Gv)
        baryv <- apply(trvertices, 2, mean)
        # euclidian dstances to Gv (dB)
        distbaryv <- rep(0, nbSpecies)
        for (j in seq_len(nbSpecies)) distbaryv[j] <- (sum((spectraits[j, ] - baryv)^2) ) ^0.5
        # mean of dB values
        meandB <- mean(distbaryv)
        # deviations to mean of db
        devdB <- distbaryv - meandB
        # computation of FDiv
        FDiv <- (sum(devdB/nbSpecies) + meandB) / (sum(abs(devdB/nbSpecies)) + meandB)
      }
    }
  }
  if ('FEve' %in% FDmetric){
    # computation of minimum spanning tree and conversion of the 'mst' matrix into 'dist' class
    tr.dist <- stats::dist(spectraits)
    linkmst <- ape::mst(tr.dist)
    mstvect <- as.dist(linkmst)
    # computation of EW for the (nbSpecies - 1) segments to link the nbSpecies points
    EW <- rep(0, nbSpecies - 1)
    flag <- 1
    for (m in seq_len((nbSpecies - 1) * nbSpecies / 2)) {
      if (mstvect[m] != 0) {
        EW[flag] <- tr.dist[m]
        flag <- flag + 1
      }
    }
    # computation of the PEW and comparison with 1 / nbSpecies - 1, finally computation of FEve
    minPEW <- rep(0, nbSpecies - 1)
    OdSmO <- 1 / (nbSpecies - 1)
    for (l in seq_len(nbSpecies - 1)) {
      minPEW[l] <- min((EW[l] / sum(EW)), OdSmO)
    }
    FEve <- ((sum(minPEW)) - OdSmO) / (1 - OdSmO)
  }
  if (!is.null(p)){p()}
  return(list('FRic' = FRic, 'FEve' = FEve, 'FDiv' = FDiv))
}
