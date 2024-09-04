#' @importFrom utils read.csv
#' @importFrom ExperimentHub createHubAccessors
#' @import SpatialExperiment
.onLoad <- function(libname, pkgname) {
  fl1 <- system.file("extdata", "metadata.csv", package=pkgname)
  fl2 <- system.file("extdata", "metadata_v2.csv", package=pkgname)
  titles <- c(read.csv(fl1, stringsAsFactors=FALSE)$Title,
              read.csv(fl2, stringsAsFactors=FALSE)$Title)
  createHubAccessors(pkgname, titles)
}