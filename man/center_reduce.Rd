% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/center_reduce.R
\name{center_reduce}
\alias{center_reduce}
\title{center and reduce data matrix based on known mean and SD}
\usage{
center_reduce(X, m, sig)
}
\arguments{
\item{X}{numeric. data matrix (each column is centered/reduced)}

\item{m}{numeric. mean of each variable in the data matrix}

\item{sig}{numeric. SD of each variable in the data matrix}
}
\value{
X numeric. Centered matrix
}
\description{
center and reduce data matrix based on known mean and SD
}
