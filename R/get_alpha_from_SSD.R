#' computes alpha diversity metrics from ssd
#'
#' @param ssd numeric. spectral species distribution
#' @param nb_pix_sunlit numeric.
#' @param alpha_metrics list. alpha diversity metrics: richness, shannon, simpson
#' @param pcelim numeric. minimum proportion of pixels to consider spectral species
#' @param hill_order numeric. Hill order
#' @param p progressor object
#'
#' @return Shannon index correspnding to the distribution
#' @export

get_alpha_from_ssd <- function(ssd, nb_pix_sunlit, alpha_metrics = 'shannon',
                               pcelim = 0.02, hill_order = 1, p = NULL){
  ClusterID <- as.numeric(names(ssd))
  if (length(ClusterID)==0)
    ClusterID <- as.numeric(colnames(ssd))

  if (pcelim > 0) {
    KeepSS <- which(ssd >= pcelim)
    if (max(ssd)>1)
      KeepSS <- which(ssd >= pcelim * nb_pix_sunlit)
    ClusterID <- ClusterID[KeepSS]
    ssd <- ssd[KeepSS]
  }
  richness <- shannon <- simpson <- fisher <- hill <- NA
  if ('hill' %in% alpha_metrics)
    hill <- get_Hill(Distrib = ssd, q = hill_order)
  if ('richness' %in% alpha_metrics)
    richness <- length(ssd)
  if ('shannon' %in% alpha_metrics)
    shannon <- get_Shannon(ssd)
  if ('simpson' %in% alpha_metrics)
    simpson <- get_Simpson(ssd)
  if (!is.null(p))
    p()
  return(list('richness' = richness, 'shannon' = shannon,
              'simpson' = simpson, 'hill' = hill))
}
