###############################################################
# Script to process Amancherla et al (2025) heart transplant rejection data
# Shreya Rao, updated Mar 2025
###############################################################

# For more details on this dataset see:
# Amancherla et al. (2025): https://doi.org/10.1101/2025.02.28.640852 

library(Seurat)
library(SpatialExperiment)
library(SingleCellExperiment)

heartSPE = readRDS("GSE290577_heart_spatial_obj.rds") |> 
  as.SingleCellExperiment()

factor_cols = c("Sample", "slide", "ct_second_pass", "patient_id", 
                "biopsy_antibody_grading", "biopsy_rejection_type", 
                "biopsy_timing", "patient_antibody_grading", "sex", 
                "donor_sex", "patient_rejection_type", "treatment", 
                "patient_biopsy", "lineage", "race", "resolution_broad")

heartSPE@colData[, factor_cols] = lapply(heartSPE@colData[, factor_cols], factor)

heartSPE = SpatialExperiment(
  assays = assays(heartSPE), 
  colData = colData(heartSPE), 
  rowData = rowData(heartSPE),
  spatialCoords = as.matrix(colData(heartSPE)[, c("adj_x_centroid", "adj_y_centroid")]))

spatialCoordsNames(heartSPE) = c("x", "y")

saveRDS(heartSPE, file = "spe_Amancherla_2025.rds")
