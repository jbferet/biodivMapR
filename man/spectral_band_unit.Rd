% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_ImageProcess.R
\name{spectral_band_unit}
\alias{spectral_band_unit}
\title{Check spectral band units and convert from nanometer to micrometer or from
micrometer to nanometer}
\usage{
spectral_band_unit(HDR, Spectral_Bands)
}
\arguments{
\item{HDR}{list. header file for a raster}

\item{Spectral_Bands}{numeric. spectral band of interest.}
}
\value{
Image numeric. full size in 3 dimensions
}
\description{
Check spectral band units and convert from nanometer to micrometer or from
micrometer to nanometer
}
