% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapperBig_PCA.R
\name{wrapperBig_PCA}
\alias{wrapperBig_PCA}
\title{this function aims at applying PCA on a raster or list of rasters
in combination with the function apply_bigRaster}
\usage{
wrapperBig_PCA(input_data, input_args)
}
\arguments{
\item{input_data}{list. Image data chunk and corresponding mask if available}

\item{input_args}{list. PCA model and associated parameters required by bigRaster}
}
\value{
list. PCA transformed spectral information
}
\description{
this function aims at applying PCA on a raster or list of rasters
in combination with the function apply_bigRaster
}
