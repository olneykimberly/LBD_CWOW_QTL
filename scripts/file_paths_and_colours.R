#----------------- Libraries
.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()
library(Matrix, lib.loc = "/usr/local/lib/R/site-library")
library(dplyr)
library(data.table)
library(edgeR)
library(GenomeInfoDb, lib.loc = "/usr/local/lib/R/site-library")
library(GenomicFeatures, lib.loc = "/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs")
library(DGEobj.utils)
library(sva)
library(MatrixEQTL)
library(tidyverse)
#--- functions 
saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}
