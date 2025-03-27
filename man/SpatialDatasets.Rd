\name{SpatialDatasets}
\docType{data}
\alias{SpatialDatasets}
\alias{spe_Keren_2018}
\alias{Ferguson_Images}
\alias{spe_Ferguson_2022}
\alias{spe_Schurch_2020}
\alias{spe_Ali_2020}
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

The datasets are from several different platforms including IMC, MIBI-TOF, CODEX, Xenium, CosMx and MERFISH. They have been sourced from various publicly available sources. 

Additional examples and documentation are provided in the package vignette.

\emph{Datasets}

The \code{SpatialDatasets} package contains the following datasets:

\itemize{
\item spe_Keren_2018 (MIBI-TOF): A study on triple negative breast cancer 
containing 40 samples measured using MIBI-TOF published by [Keren et al. (2018)](https://doi.org/10.1016/j.cell.2018.08.039). 

\item Ferguson_Images (IMC): A study on head and neck cutaneous squamous 
cell carcinoma containing 44 samples measured using IMC published by [Ferguson et al. (2022)](https://doi.org/10.1158/1078-0432.CCR-22-1332). 

\item spe_Ferguson_2022 (IMC): A study on head and neck cutaneous squamous 
cell carcinoma containing 44 samples measured using IMC published by [Ferguson et al. (2022)](https://doi.org/10.1158/1078-0432.CCR-22-1332). 

\item spe_Schurch_2020 (CODEX): A study on advanced colorectal cancer 
containing 140 samples measured using CODEX published by [Schurch et al. (2020)](https://doi.org/10.1016/j.cell.2020.07.005). 

\item spe_Ali_2020 (IMC): A study on breast cancer 
containing 483 samples measured using IMC published by [Ali et al. (2020)](https://doi.org/10.1038/s43018-020-0026-6).

\item spe_Amancherla_2025 (Xenium): A study on heart transplant rejection containing 62 (49 adult and 13 pediatric) samples measured with Xenium 
published by [Amancherla et al. (2025)](https://doi.org/10.1101/2025.02.28.640852).

\item spe_Vannan_2025 (Xenium): A study on pulmonary fibrosis containing 35 lung samples 
measured with Xenium published by [Vannan et al. (2025)](https://doi.org/10.1038/s41588-025-02080-x).
}
\value{
The datasets will be a SpatialExperiment or MoleculeExperiment.
}
\examples{
# load using dataset name
spe <- spe_Keren_2018()
spe
}