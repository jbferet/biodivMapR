#' identifies which alpha/beta metrics should be used for next processes

#' @param obs_criterion character. richness, shannon, simpson, hill, BC
#'
#' @return list including alpha_metrics and beta
#' @export

which_alpha_beta <- function(obs_criterion){
  alphamet <- c('richness', 'shannon', 'simpson', 'hill')
  betamet <- 'BC'
  # if computation of functional metrics required
  alpha_metrics <- alphamet[which(alphamet %in% obs_criterion)]
  if (length(alpha_metrics)==0)
    alpha_metrics <- NULL
  # computation of beta diversity required?
  betametrics <- betamet[which(betamet %in% obs_criterion)]
  if (length(betametrics)==0)
    getBeta <- FALSE
  if (length(betametrics)>0)
    getBeta <- TRUE

  return(list('alpha_metrics' = alpha_metrics, 'getBeta' = getBeta))
}
