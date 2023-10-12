###############################################################
# Script to process ALi et al (2020) METABRIC cohort breast cancer data
# Ellis Patrick, updated Oct 2023
###############################################################

# For more details on this dataset see:
# Ali et al. (2020): https://doi.org/10.1038/s43018-020-0026-6

# The raw data was downloaded from the Image Data Resource:
# https://idr.openmicroscopy.org/webclient/?show=project-1302
# And the supplementary data on the manuscript:
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1007-8/MediaObjects/41586_2019_1007_MOESM7_ESM.txt



# Load in data
datasetFolder <- "/dski/nobackup/biostat/datasets/spatial/IMC_BreastCancer_metabric_Ali2020"
IMC <- read.csv(file.path(datasetFolder, "Data", "single_cell_data.csv"))


# Creating patient clinical characteristics dataframe

clinical <- read.delim(file.path(datasetFolder, "Data", "NIHMS1520488-supplement-Supp_Table_5.txt"),
                       check.names = FALSE)
clinical$metabricId <- clinical$METABRIC.ID # Column name identical for joining to omics data.
clinical$timeRFS <- apply(clinical[, c("T", "TLR", "TDR")], 1, min)
clinical$eventRFS <- apply(clinical[, c("DeathBreast", "LR", "DR")], 1, max)
rownames(clinical) <- clinical[, "METABRIC.ID"]

# There is a small amount of missing data in the clinical table. 
# Impute it using random forest. 

set.seed(51773)
library(missRanger)
clinical <- missRanger(clinical, . - MATCHED.NORMAL.METABRIC.ID ~ . - METABRIC.ID - MATCHED.NORMAL.METABRIC.ID - Cohort - Date.Of.Diagnosis - Complete.Rec.History - metabricId)

clinical <- clinical |>
  select(-c("METABRIC.ID", "MATCHED.NORMAL.METABRIC.ID", "Cohort"))

# Subset to samples in common to RNA arrays, IMC and clinical data.
commonIDs <- Reduce(intersect, list(clinical$metabricId, IMC$metabricId))

clinical <- clinical[match(commonIDs, clinical$metabricId), ]
IMC <- IMC[IMC$metabricId %in% commonIDs,]

# Some patients have two or three images. Arbitrarily pick one image per 
# sample for those that have two or three images.

patientIDToImgID <- by(IMC, IMC$metabricId, function(patientColData)
{
  unique(patientColData$ImageNumber)[1]
})[commonIDs]
IMC <- IMC[IMC$ImageNumber %in% patientIDToImgID,]

# Define marker matrix
markerData <- IMC |>
  select(-c("file_id", "metabricId", "core_id", "ImageNumber", "ObjectNumber",
            "Location_Center_X", "Location_Center_Y", "SOM_nodes", "pg_cluster", 
            "description")) |>
  t() |>
  data.frame()

colnames(markerData) = seq_len(ncol(markerData))

# Define colData
columnData <- IMC |>
  select(c("file_id", "metabricId", "core_id", "ImageNumber", "ObjectNumber",
           "Location_Center_X", "Location_Center_Y", "SOM_nodes", "pg_cluster", 
           "description"))

rownames(columnData) = seq_len(nrow(columnData))


# Incorporating patient characteristics
columnData <- columnData |>
  left_join(clinical, by = c("metabricId"))

# Define spatial matrix
spatialData <- columnData |>
  select(c("Location_Center_X", "Location_Center_Y")) |>
  as.matrix()


# SingleCellExperiment
spe_Ali_2020 = SpatialExperiment(
  list(intensities = markerData),
  colData = columnData,
  spatialCoords = spatialData
)

saveRDS(spe_Ali_2020, file = "spe_Ali_2020.rds")