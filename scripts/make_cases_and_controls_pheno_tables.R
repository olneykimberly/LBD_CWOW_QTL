# Load phenotype data in R
df <- read.table("../snp_array/covariates_and_phenotype_files/covariates_and_phenotypes_imputed.txt", header=T)

# 1. Controls and PA  only (Pheno == 1)
controls <- df[df$Pheno == 1, c("FID", "IID")]
write.table(controls, "../snp_array/covariates_and_phenotype_files/controls_and_PA.txt", row.names=F, col.names=T, quote=F)

# 2. Cases which include AD and LBD (Pheno == 2)
cases <- df[df$Pheno == 2, c("FID", "IID")]
write.table(cases, "../snp_array/covariates_and_phenotype_files/cases_AD_LBD.txt", row.names=F, col.names=T, quote=F)