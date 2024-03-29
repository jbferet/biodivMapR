% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_ImageProcess.R
\name{split_image}
\alias{split_image}
\title{defines the number of pieces resulting from image split}
\usage{
split_image(HDR, LimitSizeGb = FALSE)
}
\arguments{
\item{HDR}{list. Header file}

\item{LimitSizeGb}{numeric. maximum size of individual pieces of an image (in Gb)}
}
\value{
nbPieces number of pieces
}
\description{
defines the number of pieces resulting from image split
}
