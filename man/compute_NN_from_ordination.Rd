% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_NN_from_ordination.R
\name{compute_NN_from_ordination}
\alias{compute_NN_from_ordination}
\title{compute the nearest neighbors among kernels}
\usage{
compute_NN_from_ordination(MatBC, knn, PCoA_train)
}
\arguments{
\item{MatBC}{matrix of BC dissimilarity between the kernels excluded from Ordination (rows)}

\item{knn}{numeric. number of neighbors}

\item{PCoA_train}{numeric. PCoA coordinates computed from dissimilarity matrix}
}
\value{
Ordin_est estimated NMDS position based on nearest neighbors from NMDS
}
\description{
compute the nearest neighbors among kernels
}
