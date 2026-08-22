#' computes alpha diversity metrics from ssd_mat corresponding to iterations
#'
#' @param ssd_mat numeric. spectral species distribution as a matrix:
#' species = columns, iterations = rows
#' @param alpha_metrics list. alpha diversity metrics: richness, shannon, simpson, hill
#' @param Hill_order numeric. Hill order
#'
#' @return alpha diversity metrics corresponding to the distribution
#' @importFrom vegan fisher.alpha
#' @export

get_alpha <- function(ssd_mat, alpha_metrics, Hill_order = 1){

  compute_richness <- 'richness' %in% alpha_metrics
  compute_shannon <- 'shannon' %in% alpha_metrics
  compute_simpson <- 'simpson' %in% alpha_metrics
  compute_hill <- 'hill' %in% alpha_metrics
  alpha <- compute_alpha_diversity_metrics(ssd_mat,
                                           compute_richness = compute_richness,
                                           compute_shannon = compute_shannon,
                                           compute_simpson = compute_simpson,
                                           compute_hill = compute_hill,
                                           q = Hill_order)
  return(alpha)
}
