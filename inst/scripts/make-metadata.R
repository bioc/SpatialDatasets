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
  DispatchClass = "Rds", 
  stringsAsFactors = FALSE
)


# metadata for individual datasets

df_spe_Keren_2018 <- cbind(
  df_all, 
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


# combine and save as .csv spreadsheet file

df_combined <- rbind(
  df_spe_Keren_2018
)

write.csv(df_combined, file = "metadata.csv", row.names = FALSE)

