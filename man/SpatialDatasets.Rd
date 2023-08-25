\name{SpatialDatasets}
\docType{data}
\alias{SpatialDatasets}
\alias{spe_Keren_2018}
\title{
A collection of publically available spatial omics datasets.
}
\description{
This is a collection of publically available spatial omics datasets.
  Where possible we have curated these datasets as either SpatialExperiments, 
  MoleculeExperiments or CytoImageLists and included annotations of the sample
  characteristics. 
}
\details{
The \code{SpatialDatasets} package contains a collection of spatially-resolved 
omics datasets for use in examples, demonstrations, and tutorials.

The datasets are from several different platforms including IMC, MIBI-TOF, Xenium, CosMx and MERFISH. They have been sourced from various publicly available sources. 

Additional examples and documentation are provided in the package vignette.

\emph{Datasets}

The \code{SpatialDatasets} package contains the following datasets:

\itemize{
\item spe_Keren_2018 (MIBI-TOF): A study on triple negative breast cancer 
containing 40 samples measured using MIBI-TOF published by [Keren et al. (2018)](https://doi.org/10.1016/j.cell.2018.08.039). 
}
}
\examples{
# load using dataset name
spe <- spe_Keren_2018()
spe
}