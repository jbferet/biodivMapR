% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lib_PerformPCA.R
\name{write_PCA_raster}
\alias{write_PCA_raster}
\title{writes an ENVI image corresponding to PCA}
\usage{
write_PCA_raster(
  Input_Image_File,
  Input_Mask_File,
  PCA_Path,
  PCA_model,
  Spectral,
  Nb_PCs,
  Continuum_Removal,
  TypePCA,
  nbCPU = 1,
  MaxRAM = 0.25
)
}
\arguments{
\item{Input_Image_File}{path for the raster on which PCA is applied}

\item{Input_Mask_File}{path for the corresponding mask}

\item{PCA_Path}{path for resulting PCA}

\item{PCA_model}{PCA model description}

\item{Spectral}{spectral information to be used in the image}

\item{Nb_PCs}{number of components kept in the resulting PCA raster}

\item{Continuum_Removal}{boolean. If TRUE continuum removal is performed.}

\item{TypePCA}{PCA, SPCA, NLPCA}

\item{nbCPU}{number of CPUs to process data}

\item{MaxRAM}{max RAM when initial image is read (in Gb)}
}
\value{
None
}
\description{
writes an ENVI image corresponding to PCA
}
