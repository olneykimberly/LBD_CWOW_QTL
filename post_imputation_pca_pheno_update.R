library(data.table)

setwd("/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/scripts/")

pheno <- fread("../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes.txt")
pcs <- fread("../snp_array/final_gwas_dataset/CWOW_TOPMED_final_postQC.eigenvec", header = FALSE)

setnames(pcs, c("FID", "IID", paste0("PC", 1:(ncol(pcs) - 2))))

covar <- merge(
  pcs[, .(FID, IID, PC1, PC2, PC3, PC4, PC5)],
  pheno[, .(
    FID,
    IID,
    Age,
    Sex,
    Pheno,
    Braak = Braak.NFT,
    Thal = Thal.amyloid,
    CingLB = Cing.LB
  )],
  by = c("FID", "IID"),
  all.x = TRUE
)

fwrite(
  covar,
  "../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed_final.txt",
  sep = "\t"
)