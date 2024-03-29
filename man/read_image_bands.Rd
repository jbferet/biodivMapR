% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_ImageProcess.R
\name{read_image_bands}
\alias{read_image_bands}
\title{read specific image bands from image}
\usage{
read_image_bands(ImPath, HDR, ImBand)
}
\arguments{
\item{ImPath}{Path of the image to read}

\item{HDR}{Header for the image}

\item{ImBand}{Bands to be read}
}
\value{
Image_Subset information corresponding to ImBand
}
\description{
read specific image bands from image
}
