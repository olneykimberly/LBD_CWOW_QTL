.libPaths(c("/home/kolney/R/x86_64-pc-linux-gnu-library/"))
library(data.table)
setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")

meta <- fread("../metadata/n579_metadata.txt")
# Check available group labels
table(meta$TYPE, useNA = "ifany")

# controls + PA vs AD
pheno_ctrlPA_vs_AD <- meta[
  TYPE %in% c("CONTROL", "PA", "AD"),
  .(
    FID,
    IID,
    Pheno_ctrlPA_vs_AD = fifelse(TYPE %in% c("CONTROL", "PA"), 1L, 2L)
  )
]

fwrite(
  pheno_ctrlPA_vs_AD,
  "../snp_array/covariates_and_phenotype_files/Pheno_controlsPA_vs_AD.txt",
  sep = "\t"
)

# controls + PA vs LBD
pheno_ctrlPA_vs_LBD <- meta[
  TYPE %in% c("CONTROL", "PA", "LBD"),
  .(
    FID,
    IID,
    Pheno_ctrlPA_vs_LBD = fifelse(TYPE %in% c("CONTROL", "PA"), 1L, 2L)
  )
]

fwrite(
  pheno_ctrlPA_vs_LBD,
  "../snp_array/covariates_and_phenotype_files/Pheno_controlsPA_vs_LBD.txt",
  sep = "\t"
)

# AD vs LBD (exploratory)
pheno_AD_vs_LBD <- meta[
  TYPE %in% c("AD", "LBD"),
  .(
    FID,
    IID,
    Pheno_AD_vs_LBD = fifelse(TYPE == "AD", 1L, 2L)
  )
]

fwrite(
  pheno_AD_vs_LBD,
  "../snp_array/covariates_and_phenotype_files/Pheno_AD_vs_LBD.txt",
  sep = "\t"
)

# Quick checks
cat("\ncontrols+PA vs AD:\n")
print(table(pheno_ctrlPA_vs_AD$Pheno_ctrlPA_vs_AD))

cat("\ncontrols+PA vs LBD:\n")
print(table(pheno_ctrlPA_vs_LBD$Pheno_ctrlPA_vs_LBD))

cat("\nAD vs LBD:\n")
print(table(pheno_AD_vs_LBD$Pheno_AD_vs_LBD))
