% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_ImageProcess.R
\name{extract.big_raster}
\alias{extract.big_raster}
\title{Extract bands of sparse pixels in image data cube}
\usage{
extract.big_raster(ImPath, rowcol, MaxRAM = 0.5)
}
\arguments{
\item{ImPath}{character. Path to the image cube}

\item{rowcol}{matrix or data.frame with two columns: row, col.
If columns are not named, 1st=row, 2nd=col.}

\item{MaxRAM}{numeric. Maximum memory use at block reading.
It constrains the maximum number rows of a block}
}
\value{
matrix. Rows are corresponding to the samples, columns are the bands.
}
\description{
Extract bands of sparse pixels in image data cube
}
