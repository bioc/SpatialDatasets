# ---------------------------
# Create metadata spreadsheet
# ---------------------------

# metadata for all datasets

df_all <- data.frame(
  BiocVersion = "3.18", 
  Genome = NA, 
  SourceVersion = NA, 
  Coordinate_1_based = NA, 
  DataProvider = NA, 
  Maintainer = "Ellis Patrick <ellis.patrick@sydney.edu.au>", 
  stringsAsFactors = FALSE
)


# metadata for individual datasets

# spe_Keren_2018

df_spe_Keren_2018 <- cbind(
  df_all, 
  DispatchClass = "Rds", 
  Title = "spe_Keren_2018", 
  Description =  "A study on triple negative breast cancer containing samples measured using MIBI-TOF", 
  SourceUrl = "https://www.angelolab.com/mibi-data", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataPath = "SpatialDatasets/spe_Keren_2018.rds", 
  RDataClass = "SpatialExperiment", 
  SourceType = "TIFF", 
  stringsAsFactors = FALSE
)

# Ferguson_Images

df_Ferguson_Images <- cbind(
  df_all, 
  DispatchClass = "FilePath", 
  Title = "Ferguson_Images", 
  Description =  "A study on head and neck cutaneous squamous cell carcinomas containing samples measured using IMC", 
  SourceUrl = "https://ellispatrick.github.io/", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataPath = "SpatialDatasets/Ferguson_Images", 
  RDataClass = "character", 
  SourceType = "TIFF", 
  stringsAsFactors = FALSE
)

# spe_Ferguson_2022

df_spe_Ferguson_2022 <- cbind(
  df_all, 
  DispatchClass = "rda",
  Title = "spe_Ferguson_2022", 
  Description =  "A study on head and neck cutaneous squamous cell carcinomas containing samples measured using IMC", 
  SourceUrl = "https://ellispatrick.github.io/", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataPath = "SpatialDatasets/spe_Ferguson_2022/spe_Ferguson_2022.rda", 
  RDataClass = "SpatialExperiment", 
  SourceType = "TIFF", 
  stringsAsFactors = FALSE
)

# spe_Schurch_2020

df_spe_Schurch_2020 <- cbind(
  df_all, 
  DispatchClass = "Rds", 
  Title = "spe_Schurch_2020", 
  Description =  "A study on advanced stage colorectal cancer containing samples measured using CODEX", 
  SourceUrl = "https://data.mendeley.com/public-files/datasets/mpjzbtfgfr/files/c24351b3-76d7-444f-9edf-0246356b0c78/file_downloaded", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataPath = "SpatialDatasets/spe_Schurch_2020/spe_Schurch_2020.rds", 
  RDataClass = "SpatialExperiment", 
  SourceType = "TIFF", 
  stringsAsFactors = FALSE
)

# spe_Ali_2020

df_spe_Ali_2020 <- cbind(
  df_all, 
  DispatchClass = "Rds", 
  Title = "spe_Ali_2020", 
  Description =  "A study on breast cancer containing samples measured using IMC", 
  SourceUrl = "https://idr.openmicroscopy.org/webclient/?show=project-1302", 
  Species = "Homo sapiens", 
  TaxonomyId = "9606", 
  RDataPath = "SpatialDatasets/spe_Ali_2020/spe_Ali_2020.rds", 
  RDataClass = "SpatialExperiment", 
  SourceType = "TIFF", 
  stringsAsFactors = FALSE
)

# combine and save as .csv spreadsheet file

df_combined <- rbind(
  df_spe_Keren_2018,
  df_Ferguson_Images,
  df_spe_Ferguson_2022,
  df_spe_Schurch_2020,
  df_spe_Ali_2020
)

write.csv(df_combined, file = "metadata_v2.csv", row.names = FALSE)

