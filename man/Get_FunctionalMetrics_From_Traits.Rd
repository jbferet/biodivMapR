% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_MapFunctionalDiversity.R
\name{Get_FunctionalMetrics_From_Traits}
\alias{Get_FunctionalMetrics_From_Traits}
\title{Prepare for the computation of the functional diversity metrics}
\usage{
Get_FunctionalMetrics_From_Traits(
  ReadWrite,
  Functional_File,
  Selected_Features,
  MinMaxRaster,
  HDR,
  FunctIN_Format,
  ImgFormat,
  window_size,
  MinSun,
  FDmetric = c("FRic", "FEve", "FDiv"),
  nbCPU = 1
)
}
\arguments{
\item{ReadWrite}{numeric. bytes coordinates for each read and write}

\item{Functional_File}{character. path for the raster file to get functional traits from}

\item{Selected_Features}{numeric. features to be used from Functional_File}

\item{MinMaxRaster}{numeric. min and max values to be used for data standardization}

\item{HDR}{list. header file}

\item{FunctIN_Format}{list. image format (bytes, data type, etc)}

\item{ImgFormat}{character. define if image is 2D or 3D}

\item{window_size}{numeric. window size to compute metrics from}

\item{MinSun}{numeric. minimum sun exposition to keep window as valid data (ratio of sunlit pixels / total pixels in the window)}

\item{FDmetric}{character. enumerate diversity metrics to be computed}

\item{nbCPU}{numeric. Number of CPUs to use in parallel.}
}
\value{
FDmetrics
}
\description{
Prepare for the computation of the functional diversity metrics
}
