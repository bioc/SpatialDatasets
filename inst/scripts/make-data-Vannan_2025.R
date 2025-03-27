###############################################################
# Script to process Vannan et al. (2025) pulmonary fibrosis lung data
# Shreya Rao, updated Mar 2025
###############################################################

# For more details on this dataset see:
# Vannan et al. (2025): https://doi.org/10.1038/s41588-025-02080-x

library(Seurat)
library(SpatialExperiment)
library(SingleCellExperiment)

lungSPE = readRDS("GSE250346_Seurat_GSE250346_CORRECTED_SEE_RDS_README_082024.rds") |> 
  as.SingleCellExperiment()

factor_cols = c("sample", "patient", "sample_type", "sample_affect", "disease_status",
                "tma", "run", "final_CT", "final_lineage", "CNiche", "TNiche")

lungSPE@colData[, factor_cols] = lapply(lungSPE@colData[, factor_cols], factor)

lungSPE = SpatialExperiment(
  assays = assays(lungSPE), 
  colData = colData(lungSPE), 
  rowData = rowData(lungSPE),
  spatialCoords = as.matrix(colData(lungSPE)[, c("adj_x_centroid", "adj_y_centroid")]))

spatialCoordsNames(lungSPE) = c("x", "y")

saveRDS(lungSPE, file = "spe_Vannan_2025.rds")