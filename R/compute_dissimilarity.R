#' compute dissimilarity matrix corresponding to two distinct matrices, each of
#' them including a list of samples (rows) defined by their spectral species (columns)
#' metrics available: bray, brayturn, simpson_diss, jaccard, jaccardturn, sorensen
#
#' @param A A numeric matrix of size N x p.
#' @param B A numeric matrix of size M x p.
#' @param beta_metrics list. beta diversity metrics
#
#' @return mat_diss dissimilarity matrix
#' @export

compute_dissimilarity <- function(A, B, beta_metrics) {
  mat_diss <- list()
  for (beta in beta_metrics){
    if (beta == 'bray')
      mat_diss[[beta]] <- bray_curtis_dissimilarity(A = A, B = B)
    if (beta == 'brayturn')
      mat_diss[[beta]] <- bray_curtis_turnover(A = A, B = B)
    if (beta == 'simpson_diss')
      mat_diss[[beta]] <- simpson_dissimilarity(A = A, B = B)
    if (beta == 'jaccard')
      mat_diss[[beta]] <- jaccard_dissimilarity(A = A, B = B)
    if (beta == 'jaccardturn')
      mat_diss[[beta]] <- jaccard_turnover(A = A, B = B)
    if (beta == 'sorensen')
      mat_diss[[beta]] <- sorensen_dissimilarity(A = A, B = B)
  }
  # mat_bc <- dissUtils::diss(ssd_list[[1]], ssd_list[[2]], method = 'braycurtis')
  return(mat_diss)
}
