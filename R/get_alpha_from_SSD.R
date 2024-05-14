#' computes alpha diversity metrics from SSD
#'
#' @param SSD numeric. spectral species distribution
#' @param nbPix_Sunlit numeric.
#' @param alphametrics list. alpha diversity metrics: richness, shannon, simpson
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#'
#' @return Shannon index correspnding to the distribution
#' @importFrom vegan fisher.alpha
#' @export

get_alpha_from_SSD <- function(SSD, nbPix_Sunlit, alphametrics = 'shannon', pcelim = 0.02){
  ClusterID <- as.numeric(names(SSD))
  if (pcelim > 0) {
    KeepSS <- which(SSD >= pcelim * nbPix_Sunlit)
    ClusterID <- ClusterID[KeepSS]
    SSD <- SSD[KeepSS]
  }
  richness <- shannon <- simpson <- fisher <- NA
  if (!is.na(match('richness',alphametrics))) richness <- length(SSD)
  if (!is.na(match('shannon',alphametrics))) shannon <- get_Shannon(SSD)
  if (!is.na(match('simpson',alphametrics))) simpson <- get_Simpson(SSD)
  if (!is.na(match('fisher',alphametrics)) & length(SSD)>1) fisher <- vegan::fisher.alpha(SSD)
  return(list('richness' = richness, 'shannon' = shannon,
              'simpson' = simpson, 'fisher' = fisher))
}
