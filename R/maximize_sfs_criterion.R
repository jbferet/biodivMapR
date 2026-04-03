#' define variable maximizing a criterion
#'
#' @param CorrSFS numeric
#' @param obs_criterion character. richness, shannon, simpson or BC
#' @param nbvars2select numeric
#'
#' @return SelVar numeric
#'
#' @export

maximize_sfs_criterion <- function(CorrSFS, obs_criterion, nbvars2select){
  # if a unique criterion is defined: select best
  criterion <- CorrSFS[[nbvars2select]][obs_criterion]
  if (length(obs_criterion)==1){
    SelVar <- which(criterion[[1]] == max(criterion,na.rm = T))
  } else if (length(obs_criterion)>1){
    # if several criterions are defined: rank for each criterion, and get minimum average rank
    ranks <- list()
    for (crit in obs_criterion){
      ranks[[crit]] <- sort(CorrSFS[[nbvars2select]][[crit]],
                            index.return = TRUE, decreasing = T)$ix
    }
    ranks <- as.data.frame(ranks)
    ranks_all <- rowSums(ranks)
    SelVar <- which(ranks_all == min(ranks_all,na.rm = T))
  }
  return(SelVar)
}
