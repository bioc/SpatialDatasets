###############################################################
# Script to process Keren et al (2018) triple negative breast cancer MIBI-TOF data
# Ellis Patrick, updated Aug 2023
###############################################################

# For more details on this dataset see:
# Keren et al. (2018): https://doi.org/10.1016/j.cell.2018.08.039

# The raw data was downloaded from the Angelo Lab website:
# https://www.angelolab.com/mibi-data
# And the supplementary data on the manuscript:
# https://www.cell.com/cms/10.1016/j.cell.2018.08.039/attachment/fde24f46-1466-4a73-b1c4-83672f7cf347/mmc2.xlsx



library(SpatialExperiment)
library(EBImage)
library(dplyr)
library(BiocParallel)
library(tidyr)
library(readr)
library(tibble)



nCores = 30
######
# Load data
######

patientData = read.csv("patient_class.csv", header = FALSE)
markerInfo = read.csv("cellData.csv")
patientChar = readxl::read_xlsx("mmc2.xlsx")

tiffs = list.files("images", full.names = TRUE) 
tiffnames = basename(tiffs) %>% parse_number()

raw_images = bplapply(tiffs, EBImage::readImage,
                      BPPARAM = MulticoreParam(workers = nCores))
names(raw_images) = tiffnames




convertTiff = function(image) {
  X <- matrix(1:nrow(image), nrow = nrow(image), ncol = ncol(image), byrow = FALSE)
  Y <- matrix(1:ncol(image), nrow = nrow(image), ncol = ncol(image), byrow = TRUE)
  
  x <- tapply(X, image, mean)
  y <- tapply(Y, image, mean)
  
  df = data.frame(x, y)
  
  df$CellID = seq_len(nrow(df)) - 1
  return(df)
}


tiffdfs = bplapply(raw_images, 
                   convertTiff,
                   BPPARAM = MulticoreParam(workers = nCores))



for(i in 1:length(tiffdfs)) {
  image = tiffdfs[[i]]
  image$imageID = as.numeric(names(tiffdfs)[[i]])
  image = remove_rownames(image)
  tiffdfs[[i]] = image
}

spatialData = tiffdfs |> 
  bind_rows() |> 
  inner_join(markerInfo,
             by = c("imageID" = "SampleID", "CellID" = "cellLabelInImage")
  ) %>%
  arrange(imageID, CellID)


tumour <- c("Keratin_Tumour", "Tumour")
bcells <- c("B_cell")
tcells <- c("dn_T_cell", "CD4_T_cell", "CD8_T_cell", "Tregs")
myeloid <- c("Dc_or_Mono", "DC", "Mono_or_Neu", "Macrophages", "Other_Immune", "Neutrophils")


# Labelling cell types
spatialData = spatialData %>%
  mutate(
    cellType = case_when(
      Group == 1 ~ "Unidentified",
      immuneGroup == 1 ~ "Tregs",
      immuneGroup == 2 ~ "CD4_T_cell",
      immuneGroup == 3 ~ "CD8_T_cell",
      immuneGroup == 4 ~ "dn_T_CD3",
      immuneGroup == 5 ~ "NK",
      immuneGroup == 6 ~ "B_cell",
      immuneGroup == 7 ~ "Neutrophils",
      immuneGroup == 8 ~ "Macrophages",
      immuneGroup == 9 ~ "DC",
      immuneGroup == 10 ~ "DC_or_Mono",
      immuneGroup == 11 ~ "Mono_or_Neu",
      immuneGroup == 12 ~ "Other_Immune",
      Group == 3 ~ "Endothelial",
      Group == 4 ~ "Mesenchymal",
      Group == 5 ~ "Tumour",
      Group == 6 ~ "Keratin_Tumour"
    )
  )



# Creating patient labels: Factor 0 is mixed tumours, 1 = compartmentalised, 2 = cold

#Patients with less than 250 immune cells (N = 6) were defined as cold. Patients with a mixing score < 0.22 (N = 15) were defined as compartmentalized and the rest of the patients (N = 20) were defined as mixed
# patientData$V2 %>% 
#   table() 


patientData = patientData %>% 
  rename("patient" = V1,
         "tumour_type" = V2) %>%
  mutate(tumour_type = case_when(tumour_type == 0 ~ "mixed",
                                 tumour_type == 1 ~ "compartmentalised",
                                 tumour_type == 2 ~ "cold") %>% 
           factor())


# Merging patient data
spatialData = left_join(spatialData, patientData, by = c("imageID" = "patient"))


# Define colData
columnData = spatialData |> 
  select("x", "y", "CellID", "imageID", "cellType", "tumour_type", "cellSize",
         "C", "tumorYN", "tumorCluster", "Group", "immuneCluster", "immuneGroup")

rownames(columnData) = seq_len(nrow(columnData))



# Incorporating patient characteristics
columnData = columnData |> 
  left_join(patientChar, by = c("imageID" = "InternalId"))


# Define marker matrix
markerData = spatialData |> 
  select(-c("x", "y", "CellID", "imageID", "cellType", "tumour_type", "cellSize",
            "C", "tumorYN", "tumorCluster", "Group", "immuneCluster", "immuneGroup")) |> 
  t() |> 
  data.frame()

colnames(markerData) = seq_len(ncol(markerData))


# SingleCellExperiment
spe_Keren_2018 = SpatialExperiment(
  list(intensities = markerData),
  colData = columnData
)


saveRDS(spe_Keren_2018, file = "spe_Keren_2018.rds")