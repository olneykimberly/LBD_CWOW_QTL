#----------------- Libraries
.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()
library(Matrix, lib.loc = "/usr/local/lib/R/site-library")
library(dplyr)
library(data.table)
library(edgeR)
library(GenomicFeatures)
library(DGEobj.utils)
library(sva)
library(MatrixEQTL)
#--- functions 
saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts")
# Read in metadata filtered metadata file 
metadata <- read.delim("../metadata/n580_metadata.txt")

# Read the raw genotype file
genotype_raw <- fread("../snp_array/TOPMED_imput/final_filtered_data.genotype.raw", sep = " ")
rownames(genotype_raw) <- genotype_raw$IID

# Remove the first six columns (FID, IID, PAT, MAT, SEX, PHENOTYPE)
library(data.table)
genotype_data <- genotype_raw[, -(1:6), with = FALSE]
rownames(genotype_data) <- genotype_raw$IID
# Transpose the data to get SNPs as rows and samples as columns
genotype_data <- t(genotype_data)
colnames(genotype_data) <- genotype_raw$IID

# Save the formatted data to a file
write.table(genotype_data, file = "../snp_array/final_filtered_genotype_data.txt", sep = " ", quote = FALSE, col.names = NA)