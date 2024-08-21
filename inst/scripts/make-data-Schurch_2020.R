###############################################################
# Script to process Schurch et al (2020) CODEX colorectcal cancer data
# Alex Qin, updated Oct 2023
###############################################################

# For more details on this dataset see:
# Schurch et al. (2020): https://doi.org/10.1016/j.cell.2020.07.005

# The raw data was downloaded from the Image Data Resource:
# https://data.mendeley.com/public-files/datasets/mpjzbtfgfr/files/c24351b3-76d7-444f-9edf-0246356b0c78/file_downloaded
# And the supplementary data on the manuscript:
# https://www.cell.com/cms/10.1016/j.cell.2020.07.005/attachment/bd873aa1-79c6-4de4-8b06-905f37bf478f/mmc2.xlsx



# Load in data
codexData <- readr::read_csv("/dskh/nobackup/biostat/datasets/spatial/CODEX_Colon_Schurch2020/Data/CRC_clusters_neighborhoods_markers.csv") %>%
  dplyr::select(-1)
lev <- unique(codexData$ClusterName)
codexData$cellType <- factor(codexData$ClusterName, levels = lev, labels = janitor::make_clean_names(lev))

# Creating marker matrix

markerData <- codexData[, grepl("-", names(codexData))]

colnames(markerData)[str_detect(colnames(markerData)," - ")] = colnames(markerData)[str_detect(colnames(markerData)," - ")] %>%
  stringr::str_split(" - ") %>% lapply(function(x) x[[1]])        

markerData = markerData %>%   
  janitor::clean_names() %>%
  as.data.frame()

rownames(markerData) <- paste0("cellIDcell_", rownames(markerData))

markerData = t(markerData)

# Define colData
clinicalCodex <- readxl::read_xlsx("/dskh/nobackup/biostat/datasets/spatial/CODEX_Colon_Schurch2020/Data/CRC_TMAs_patient_annotations.xlsx") %>% janitor::clean_names()

phenoData <- clinicalCodex %>% 
  mutate(patients = patient) %>%
  as.data.frame()

imageMatch <- codexData %>% 
  mutate(cellID = rownames(codexData)) %>%
  select(imageID = "File Name", 
         cellID,
         x = "X:X",
         y = "Y:Y",
         cellType,
         patients)


phenoData <- left_join(imageMatch, phenoData, by = "patients")  %>% 
  as.data.frame() %>%
  mutate(subject = gsub("_A|_B", "", imageID), patient = factor(patients), group = factor(group))


tissueType = phenoData %>% 
  select(imageID, cellID, la_4, la_6, diffuse_5, diffuse_7, subject = patient) %>%
  mutate(regionNumber = as.numeric(stringr::str_sub(imageID, 6, 6))) %>% 
  mutate(subSpot = ifelse(str_detect(imageID, "A"), 1, 2)) %>% 
  group_by(subject) %>% 
  mutate(regionOrder = ifelse(regionNumber == min(regionNumber), 1, 2)) %>% 
  ungroup() %>% 
  select(-regionNumber) %>% 
  mutate(tissueType = case_when(
    la_4 == 1 &  regionOrder == 1 & subSpot == 1 ~ 0, 
    la_4 == 1 &  regionOrder == 1 & subSpot == 2 ~ 1,
    la_4 == 2 &  regionOrder == 1 ~ 0,
    diffuse_5 == 2 &  regionOrder == 1 ~ 1,
    la_6 == 1 &  regionOrder == 2 & subSpot == 1 ~ 0, 
    la_6 == 1 &  regionOrder == 2 & subSpot == 2 ~ 1,
    la_6 == 2 &  regionOrder == 2 ~ 0,
    diffuse_7 == 2 &  regionOrder == 2 ~ 1
  )) %>% 
  select(imageID, tissueType, cellID)

phenoData = phenoData %>% 
  left_join(tissueType, by = c("cellID", "imageID"))

rownames(phenoData) = paste0("cellIDcell_",phenoData$cellID)

# Define spatial matrix
spatialData <- phenoData |>
  select(c("x", "y")) |>
  as.matrix()


# SingleCellExperiment
spe_Schurch_2020 <- SpatialExperiment::SpatialExperiment(
  list(intensities = markerData),
  colData = phenoData,
  spatialCoords = spatialData
)

saveRDS(spe_Schurch_2020, file = "spe_Schurch_2020.rds")