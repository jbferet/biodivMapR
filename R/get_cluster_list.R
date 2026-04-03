#' get cluster list for biodivMapR_opt_clusters

#' @param nb_clusters numeric. number of clusters
#'
#' @return list including alpha_metrics and beta
#' @export

get_cluster_list <- function(nb_clusters){
  if (length(nb_clusters)==1)
    nbClust_list <- seq(2,nb_clusters)
  if (length(nb_clusters)>1)
    nbClust_list <- nb_clusters
  return(nbClust_list)
}
