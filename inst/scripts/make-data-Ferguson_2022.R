###############################################################
# Script to process Ferguson et al (2022) head and neck cancer data
# Alex Qin, updated Oct 2023
###############################################################

# For more details on this dataset see:
# Ferguson et al. (2022): https://doi.org/10.1158/1078-0432.CCR-22-1332

library(cytomapper)
library(simpleSeg)
library(SpatialExperiment)

# Read in data
images <- cytomapper::loadImages(
  "Ferguson_Images",
  single_channel = TRUE,
  on_disk = TRUE,
  h5FilesPath = HDF5Array::getHDF5DumpDir()
)

mcols(images) <- S4Vectors::DataFrame(imageID = names(images))

cn <- channelNames(images) # Read in channel names
head(cn)

cn <- sub(".*_", "", cn) # Remove preceding letters
cn <- sub(".ome", "", cn) # Remove the .ome
head(cn)

channelNames(images) <- cn # Reassign channel names

head(names(images))

nam <- sapply(strsplit(names(images), "_"), `[`, 3)
head(nam)

names(images) <- nam # Reassigning image names
mcols(images)[["imageID"]] <- nam # Reassigning image names

masks <- simpleSeg(images,
                   nucleus = c("HH3"),
                   pca = TRUE,
                   cellBody = "dilate",
                   discSize = 3,
                   sizeSelection = 20,
                   transform = "sqrt",
                   tissue = c("panCK", "CD45", "HH3")
)

# Summarise the expression of each marker in each cell
cells <- cytomapper::measureObjects(masks,
                                    images,
                                    img_id = "imageID",
                                    return_as = "spe",
                                    BPPARAM = BPPARAM)

spatialCoordsNames(cells) <- c("x", "y")

clinicalData <- read.csv("ferguson_clinical.csv")

rownames(clinicalData) <- clinicalData$imageID
clinicalData <- clinicalData[names(images), ]

colData(cells) <- cbind(colData(cells), clinicalData[cells$imageID, ])

