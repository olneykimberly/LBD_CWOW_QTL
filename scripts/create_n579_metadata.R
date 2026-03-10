.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))
library(data.table)

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")

metadata <- fread("../metadata/n580_metadata.txt", sep = "\t")
fam <- fread("../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC.fam", header = FALSE)
setnames(fam, c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

# keep only metadata rows present in final PLINK file
metadata_n579 <- metadata[IID %in% fam$IID]

# reorder to match PLINK sample order
metadata_n579 <- metadata_n579[match(fam$IID, metadata_n579$IID)]

# sanity checks
stopifnot(nrow(metadata_n579) == nrow(fam))
stopifnot(all(metadata_n579$IID == fam$IID))

fwrite(metadata_n579, "../metadata/n579_metadata.txt", sep = "\t")